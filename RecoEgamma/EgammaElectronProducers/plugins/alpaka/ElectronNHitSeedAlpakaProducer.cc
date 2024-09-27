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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/ParametrizedEngine/interface/ParabolicParametrizedMagneticField.h"

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

// testing 
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/PerpendicularBoundPlaneBuilder.h"


#include "TrackingTools/TrajectoryState/interface/ftsFromVertexToPoint.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ftsFromVertexToPointPortable.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/alpaka/helixPropagator.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/alpaka/helixBarrelPlaneCrossingByCircle.h"

#include "DataFormats/EgammaReco/interface/Plane.h"

using Vector3f = Eigen::Matrix<float, 3, 1>;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

	class ElectronNHitSeedAlpakaProducer : public global::EDProducer<> {
	public:
		ElectronNHitSeedAlpakaProducer(const edm::ParameterSet& pset): 
			deviceToken_{produces()},
			size_{pset.getParameter<int32_t>("size")},
			initialSeedsToken_(consumes(pset.getParameter<edm::InputTag>("initialSeeds"))),
			beamSpotToken_(consumes(pset.getParameter<edm::InputTag>("beamSpot"))),
			magFieldToken_(esConsumes())
			{
				superClustersTokens_ = consumes(pset.getParameter<edm::InputTag>("superClusters"));
			}

		void produce(edm::StreamID sid, device::Event& event, device::EventSetup const& iSetup) const override {

			auto vprim_ = event.get(beamSpotToken_).position();
			std::cout<<"vprim "<<vprim_<<std::endl;

			// Get MagField ESProduct for comparing 
			auto const& magField = iSetup.getData(magFieldToken_);


			PropagatorWithMaterial backwardPropagator_ = PropagatorWithMaterial(oppositeToMomentum, 0.000511, &magField);


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
			auto& viewSeeds = hostProductSeeds.view();

			////////////////////////////////////////////////////////////////////
			// Fill in SOAs
			// Technically should separate in different producers that create the SoAs

			i = 0;
			printf("Printed from host : \n");
	        for (auto& superClusRef : event.get(superClustersTokens_)) 
			{
				viewSCs[i].scSeedTheta() =  superClusRef->seed()->position().theta();
				viewSCs[i].scPhi() = superClusRef->position().phi();
				viewSCs[i].scR() = superClusRef->position().r();
				viewSCs[i].scEnergy() = superClusRef->energy();
				++i;
			}

			i = 0;
			for (auto& initialSeedRef : event.get(initialSeedsToken_)) 
			{	
				viewSeeds[i].nHits() = initialSeedRef.nHits();
				for (size_t iHit = 0;iHit < initialSeedRef.nHits();iHit++) {		
					auto const& recHit = *(initialSeedRef.recHits().begin() + iHit);  
					//viewSeeds[[i][iHit]].detectorID() = recHit.geographicalId().subdetId() == PixelSubdetector::PixelBarrel ? 1: 0;
					//viewSeeds[i][iHit].isValid() = recHit.isValid();
				}
				++i;
			}

			////////////////////////////////////////////////////////////////////


	        for (auto& superClusRef : event.get(superClustersTokens_)) 
			{
				float x = superClusRef->position().r() * sin(superClusRef->seed()->position().theta()) * cos(superClusRef->position().phi());
				float y = superClusRef->position().r() * sin(superClusRef->seed()->position().theta()) * sin(superClusRef->position().phi());
				float z = superClusRef->position().r() * cos(superClusRef->seed()->position().theta());
				GlobalPoint center(x, y, z);
				float theMagField = magField.inTesla(center).mag();
				Vector3f position{x,y,z};
  				GlobalPoint vprim(vprim_.x(), vprim_.y(), vprim_.z());
				Vector3f vertex{vprim.x(),vprim.y(),vprim.z()};
				GlobalPoint sc(GlobalPoint::Polar(superClusRef->seed()->position().theta(),  //seed theta
                                                superClusRef->position().phi(),    //supercluster phi
                                                superClusRef->position().r()));    //supercluster r
			
				auto freeTS = trackingTools::ftsFromVertexToPoint(magField, sc, vprim, superClusRef->energy(), 1);
				auto initialTrajState = TrajectoryStateOnSurface(freeTS, *PerpendicularBoundPlaneBuilder{}(freeTS.position(), freeTS.momentum()));

				auto newfreeTS = ftsFromVertexToPointPortable::ftsFromVertexToPoint(position, vertex, superClusRef->energy(),1,magneticFieldParabolicPortable::magneticFieldAtPoint(position));			
				Vector3f testposition = {newfreeTS.position(0),newfreeTS.position(1),newfreeTS.position(2)};
				Vector3f testmomentum = {newfreeTS.momentum(0),newfreeTS.momentum(1),newfreeTS.momentum(2)};

				std::cout<<"  Supercluster "<< superClusRef->position().x()<<" "<< superClusRef->position().y()<<std::endl;
				std::cout<<"  Supercluster "<< x <<" "<< y <<" "<<z <<std::endl;


				auto transverseCurvature = [](const Vector3D& p, int charge, const float& magneticFieldZ) -> float {
   			 			return -2.99792458e-3f * (charge / sqrt(p(0) * p(0) + p(1) * p(1))) * magneticFieldZ;  
				};

				for (auto& initialSeedRef : event.get(initialSeedsToken_)) 
				{			


					//std::cout<<" ------------Start loop"<<std::endl;

					auto const& recHit = *(initialSeedRef.recHits().begin() + 0);  
					if(!recHit.isValid())      
						continue;

					std::cout<<" recHit.det()->surface().position() "<<recHit.det()->surface().position()<<" rotation? "<<recHit.det()->surface().rotation().z()<<std::endl;

					//std::cout<<" Print out legacy fts position "<< freeTS.position() <<"  and new implementation one : "<<testposition(0) <<" "<<testposition(1)<<" "<<testposition(2)<<std::endl;
					//std::cout<<" Print out legacy fts momentum "<< freeTS.momentum() <<"  and new implementation one : "<<testmomentum(0) <<" "<<testmomentum(1)<<" "<<testmomentum(2) <<std::endl;
					//std::cout<<" initialTrajState pos "<< initialTrajState.globalPosition()<<"  and momentum  "<<initialTrajState.globalMomentum()<<std::endl;

					double rho = 0.;
					double s = 0.;
				
					bool theSolExists = false;
					Vector3f surfPosition{recHit.det()->surface().position().x(),recHit.det()->surface().position().y(),recHit.det()->surface().position().z()};
					Vector3f surfRotation{recHit.det()->surface().rotation().z().x(),recHit.det()->surface().rotation().z().y(),recHit.det()->surface().rotation().z().z()};

					Vector3f x2{0,0,0};
					Vector3f p2{0,0,0};

					rho = transverseCurvature(testmomentum,1,magneticFieldParabolicPortable::magneticFieldAtPoint(position));
					//std::cout<<" New "<< rho <<"  and old "<<freeTS.transverseCurvature() <<std::endl;


					Propagators::helixBarrelPlaneCrossing(testposition,testmomentum,rho,oppositeToMomentum,surfPosition,surfRotation,theSolExists,x2,p2,s);
					if(!theSolExists)
						continue;

					auto state = backwardPropagator_.propagate(initialTrajState, recHit.det()->surface());
					if(!state.isValid())
						continue;


					std::cout<<" Initial: "<< state.globalParameters().position()<<"   and new "<< x2(0) <<" "<< x2(1) << " "<< x2(2)<<std::endl;

				}


				//printf("Print out legacy fts position %f and new fts position  and %d \n",freeTS.position().x(), testcharge);
				//printf("For SC i=%d Energy is :%f , theta is :%f,  r is : %f \n",i,superClusRef->energy(),superClusRef->seed()->position().theta(),superClusRef->position().r()) ;
				//printf(" view %lf ", viewSCs[i].scR());
				//printf("Magnetic field full  = %f and the parabolic approximation %f ", theMagField, magneticFieldParabolicPortable::magneticFieldAtPoint(position));
			}


			alpaka::memcpy(event.queue(), deviceProductSCs.buffer(), hostProductSCs.buffer());

			// Print the SoA 
			//algo_.print(event.queue(), deviceProductSCs);
			//algo_.matchSeeds(event.queue(), deviceProductSCs);
			event.emplace(deviceToken_, std::move(deviceProductSCs));
		}

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
			edm::ParameterSetDescription desc;
			desc.add<int32_t>("size");
			desc.add<edm::InputTag>("initialSeeds", {"hltElePixelSeedsCombined"});
			desc.add<edm::InputTag>("beamSpot", {"hltOnlineBeamSpot"});
  			desc.add<edm::InputTag>("superClusters", {"hltEgammaSuperClustersToPixelMatch"});
			descriptions.addWithDefaultLabel(desc);
		}

	private:
		const device::EDPutToken<reco::SuperclusterDeviceCollection> deviceToken_;
		const int32_t size_;
		const edm::EDGetTokenT<TrajectorySeedCollection> initialSeedsToken_;
		const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
		edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
		edm::EDGetTokenT<std::vector<reco::SuperClusterRef>> superClustersTokens_;
		SuperclusterAlgo const algo_{};
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(ElectronNHitSeedAlpakaProducer);
