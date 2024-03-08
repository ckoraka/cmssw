#include <Eigen/Core>

#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterSoA.h"
#include "DataFormats/EgammaReco/interface/EleSeedSoA.h"
#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/SuperclusterHostCollection.h"
#include "DataFormats/EgammaReco/interface/alpaka/EleSeedDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/EleSeedHostCollection.h"
#include "DataFormats/EgammaReco/interface/EleSeedHostCollection.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaServices/interface/alpaka/AlpakaService.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"

#include "SuperclusterAlgo.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE {

	class ElectronNHitSeedAlpakaProducer : public global::EDProducer<> {
	public:
		ElectronNHitSeedAlpakaProducer(const edm::ParameterSet& pset): 
			deviceToken_{produces()},
			size_{pset.getParameter<int32_t>("size")},
		    initialSeedsToken_(consumes(pset.getParameter<edm::InputTag>("initialSeeds"))),
			magFieldToken_(esConsumes())
			{
				superClustersTokens_ = consumes(pset.getParameter<edm::InputTag>("superClusters"));
			}

		void produce(edm::StreamID sid, device::Event& event, device::EventSetup const& iSetup) const override {

			// Get MagField ESProduct for comparing 
			auto const& magField = iSetup.getData(magFieldToken_);

			int i=0;
	        for (auto& superClusRef : event.get(superClustersTokens_)) {
				++i;
			}

			reco::SuperclusterHostCollection hostProductSCs{i, event.queue()};
			reco::SuperclusterDeviceCollection deviceProductSCs{i, event.queue()};

			i = 0; // To fix these make they no sense
			for (auto& initialSeedRef : event.get(initialSeedsToken_)) 
				++i;

			reco::EleSeedHostCollection hostProductSeeds{i, event.queue()};
			reco::EleSeedDeviceCollection deviceProductSeeds{i, event.queue()};

			auto& viewSCs = hostProductSCs.view();

			i = 0;
			printf("Printed from host : \n");
	        for (auto& superClusRef : event.get(superClustersTokens_)) {
				printf("For SC i=%d Energy is :%f , theta is :%f,  r is : %f \n",i,superClusRef->energy(),superClusRef->seed()->position().theta(),superClusRef->position().r()) ;
				viewSCs[i].scSeedTheta() =  superClusRef->seed()->position().theta();
				viewSCs[i].scPhi() = superClusRef->position().phi();
				viewSCs[i].scR() = superClusRef->position().r();
				printf(" view %lf ", viewSCs[i].scR());
				viewSCs[i].scEnergy() = superClusRef->energy();
				i++;
				float x = superClusRef->position().r() * sin(superClusRef->seed()->position().theta()) * cos(superClusRef->position().phi());
				float y = superClusRef->position().r() * sin(superClusRef->seed()->position().theta()) * sin(superClusRef->position().phi());
				float z = superClusRef->position().r() * cos(superClusRef->seed()->position().theta());
				GlobalPoint center(x, y, z);
				float theMagField = magField.inTesla(center).mag();
				std::cout << "Magnetic field full  = " << theMagField << std::endl;
			}
						
			alpaka::memcpy(event.queue(), deviceProductSCs.buffer(), hostProductSCs.buffer());

			// Print the SoA 
			algo_.print(event.queue(), deviceProductSCs);
			algo_.matchSeeds(event.queue(), deviceProductSCs);
			event.emplace(deviceToken_, std::move(deviceProductSCs));
		}

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
			edm::ParameterSetDescription desc;
			desc.add<int32_t>("size");
			desc.add<edm::InputTag>("initialSeeds", {"hltElePixelSeedsCombined"});
  			desc.add<edm::InputTag>("superClusters", {"hltEgammaSuperClustersToPixelMatch"});
			descriptions.addWithDefaultLabel(desc);
		}

	private:
		const device::EDPutToken<reco::SuperclusterDeviceCollection> deviceToken_;
		const int32_t size_;
		const edm::EDGetTokenT<TrajectorySeedCollection> initialSeedsToken_;
		edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
		edm::EDGetTokenT<std::vector<reco::SuperClusterRef>> superClustersTokens_;
		SuperclusterAlgo const algo_{};
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(ElectronNHitSeedAlpakaProducer);
