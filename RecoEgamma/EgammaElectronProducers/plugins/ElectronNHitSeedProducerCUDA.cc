//******************************************************************************
//
// Part of the refactorisation of of the E/gamma pixel matching for 2017 pixels
// This refactorisation converts the monolithic  approach to a series of
// independent producer modules, with each modules performing  a specific
// job as recommended by the 2017 tracker framework
//
//
// The module produces the ElectronSeeds, similarly to ElectronSeedProducer
// although with a varible number of required hits
//
//
// Author : Sam Harper (RAL), 2017
//
//*******************************************************************************

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/TrajSeedMatcher.h"


#include <cuda_runtime.h>

#include "HeterogeneousCore/CUDAUtilities/interface/host_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/device_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/host_noncached_unique_ptr.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "RecoEgamma/EgammaElectronProducers/interface/TestCUDA.h"
#include "DataFormats/EgammaReco/interface/EcalSCHeterogeneous.h"

class ElectronNHitSeedProducerCUDA : public edm::global::EDProducer<> {
public:
  explicit ElectronNHitSeedProducerCUDA(const edm::ParameterSet&);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const final;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const TrajSeedMatcher::Configuration matcherConfiguration_;

  std::vector<edm::EDGetTokenT<std::vector<reco::SuperClusterRef>>> superClustersTokens_;
  const edm::EDGetTokenT<TrajectorySeedCollection> initialSeedsToken_;
  const edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::EDGetTokenT<MeasurementTrackerEvent> measTkEvtToken_;
  const edm::EDPutTokenT<reco::ElectronSeedCollection> putToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;
};

namespace {
  int getLayerOrDiskNr(DetId detId, const TrackerTopology& trackerTopo) {
    if (detId.subdetId() == PixelSubdetector::PixelBarrel) {
      return trackerTopo.pxbLayer(detId);
    } else if (detId.subdetId() == PixelSubdetector::PixelEndcap) {
      return trackerTopo.pxfDisk(detId);
    } else
      return -1;
  }

  reco::ElectronSeed::PMVars makeSeedPixelVar(const TrajSeedMatcher::MatchInfo& matchInfo,
                                              const TrackerTopology& trackerTopo) {
    int layerOrDisk = getLayerOrDiskNr(matchInfo.detId, trackerTopo);
    reco::ElectronSeed::PMVars pmVars;
    pmVars.setDet(matchInfo.detId, layerOrDisk);
    pmVars.setDPhi(matchInfo.dPhiPos, matchInfo.dPhiNeg);
    pmVars.setDRZ(matchInfo.dRZPos, matchInfo.dRZNeg);

    return pmVars;
  }

  // GPU related vars 
  //ElectronSeedCUDA::ElectronSeedParameters eParams;

}  // namespace

ElectronNHitSeedProducerCUDA::ElectronNHitSeedProducerCUDA(const edm::ParameterSet& pset)
    : matcherConfiguration_(pset.getParameter<edm::ParameterSet>("matcherConfig"), consumesCollector()),
      initialSeedsToken_(consumes(pset.getParameter<edm::InputTag>("initialSeeds"))),
      verticesToken_(consumes(pset.getParameter<edm::InputTag>("vertices"))),
      beamSpotToken_(consumes(pset.getParameter<edm::InputTag>("beamSpot"))),
      measTkEvtToken_(consumes(pset.getParameter<edm::InputTag>("measTkEvt"))),
      putToken_{produces<reco::ElectronSeedCollection>()},
      trackerTopologyToken_{esConsumes()} {
  for (const auto& scTag : pset.getParameter<std::vector<edm::InputTag>>("superClusters")) {
    superClustersTokens_.emplace_back(consumes(scTag));
  }
}

void ElectronNHitSeedProducerCUDA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("initialSeeds", {"hltElePixelSeedsCombined"});
  desc.add<edm::InputTag>("vertices", {});
  desc.add<edm::InputTag>("beamSpot", {"hltOnlineBeamSpot"});
  desc.add<edm::InputTag>("measTkEvt", {"hltSiStripClusters"});
  desc.add<std::vector<edm::InputTag>>("superClusters", {{"hltEgammaSuperClustersToPixelMatch"}});
  desc.add<edm::ParameterSetDescription>("matcherConfig", TrajSeedMatcher::makePSetDescription());

  descriptions.add("ElectronNHitSeedProducerCUDA", desc);
}

void ElectronNHitSeedProducerCUDA::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  auto const& trackerTopology = iSetup.getData(trackerTopologyToken_);

  reco::ElectronSeedCollection eleSeeds{};
  auto CPUbeta = cms::cuda::make_host_unique<double[]>(1, cudaStreamDefault);
  //EgammaCUDA::hello_world_gpu_Wrapper();
 
  // --------- SoA DataFormat building --------------- //

  EcalSCHeterogeneous CPUecalSCs(cms::cuda::make_host_unique<EcalSC::EcalSCSoA>(cudaStreamDefault));
  EcalSCHeterogeneous GPUecalSCs(cms::cuda::make_device_unique<EcalSC::EcalSCSoA>(cudaStreamDefault));

  auto* CPUecalSCsObject = CPUecalSCs.get();
  auto* GPUecalSCsObject = GPUecalSCs.get();

  u_int16_t n_sc=0;
  for (const auto& superClustersToken : superClustersTokens_) {
    for (auto& superClusRef : iEvent.get(superClustersToken)) {
      CPUecalSCsObject->scTheta(n_sc) = superClusRef->seed()->position().theta();
      CPUecalSCsObject->scEta(n_sc) = superClusRef->seed()->position().eta();
      CPUecalSCsObject->scEnergy(n_sc) = superClusRef->seed()->position().r();
      CPUecalSCsObject->scR(n_sc) = superClusRef->energy();     
      n_sc=n_sc+1;
    }
  }

  std::cout<<" Number of SCs is : "<< n_sc <<std::endl;
  std::cout<<" Finish from Host : "<<std::endl;
  std::cout<<"sizeof(EcalSC::EcalSCSoA) : "<<sizeof(EcalSC::EcalSCSoA)<<std::endl;

  // --------- Copy list of SoA SCs on the device -----//
  cudaCheck(cudaMemcpy(GPUecalSCsObject, CPUecalSCsObject, sizeof(EcalSC::EcalSCSoA), cudaMemcpyHostToDevice));
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
      fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
  }
  EgammaCUDA::ecalScWrapper(n_sc,GPUecalSCsObject,cudaStreamDefault);
  //cudaCheck(cudaFree(GPUecalSCsObject));
  // ------------------------------------------------- //


  TrajSeedMatcher matcher{iEvent.get(initialSeedsToken_),
                          iEvent.get(beamSpotToken_).position(),
                          matcherConfiguration_,
                          iSetup,
                          iEvent.get(measTkEvtToken_)};

  int counter = 0;
  for (auto& superClusRef : iEvent.get(initialSeedsToken_)) {
    ++counter;
  }  
  std::cout<<" Iniitial seeds : "<<counter <<std::endl;

  // Loop over all super-cluster collections (typically barrel and forward are supplied separately)
  for (const auto& superClustersToken : superClustersTokens_) {
    for (auto& superClusRef : iEvent.get(superClustersToken)) {

      GlobalPoint caloPosition(GlobalPoint::Polar(superClusRef->seed()->position().theta(),  //seed theta
                                                  superClusRef->position().phi(),            //supercluster phi
                                                  superClusRef->position().r()));            //supercluster r
      int number = 0;
      for (auto const& matchedSeed : matcher(caloPosition, superClusRef->energy())) {

        ++number;
        reco::ElectronSeed eleSeed(matchedSeed.seed);
        reco::ElectronSeed::CaloClusterRef caloClusRef(superClusRef);
        eleSeed.setCaloCluster(caloClusRef);
        eleSeed.setNrLayersAlongTraj(matchedSeed.nrValidLayers);
        for (auto const& matchInfo : matchedSeed.matchInfos) {
          eleSeed.addHitInfo(makeSeedPixelVar(matchInfo, trackerTopology));
        }
        eleSeeds.emplace_back(eleSeed);
      }
      std::cout<<"Number of seeds per SC "<< number <<std::endl;
    }
  }
  
  iEvent.emplace(putToken_, std::move(eleSeeds));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronNHitSeedProducerCUDA);
