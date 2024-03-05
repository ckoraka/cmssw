#include <Eigen/Core>

#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterSoA.h"
#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/SuperclusterHostCollection.h"

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

	class SuperClusterSoAProducer : public global::EDProducer<> {
	public:
		SuperClusterSoAProducer(edm::ParameterSet const& config)
			: deviceToken_{produces()},
			size_{config.getParameter<int32_t>("size")},
			magFieldToken_(esConsumes())
			{
			for (const auto& scTag : config.getParameter<std::vector<edm::InputTag>>("getsuperclus")) {
				superClustersTokens_.emplace_back(consumes(scTag));
			}
		}

		void produce(edm::StreamID sid, device::Event& event, device::EventSetup const& iSetup) const override {

			int i=0;
			printf("Printed from host : \n");
			for (const auto& superClustersToken : superClustersTokens_) {
				for (auto& superClusRef : event.get(superClustersToken)) {
					i++;
				}
			}	

			reco::SuperclusterHostCollection hostProduct{i, event.queue()};
			reco::SuperclusterDeviceCollection deviceProduct{i, event.queue()};

			auto& view = hostProduct.view();

			// Get MagField ES product :
			auto const& magField = iSetup.getData(magFieldToken_);
			GlobalPoint center(1.0, 1.0, 1.0);
  			float theMagField = magField.inTesla(center).mag();
  			std::cout << "theMagField = " << theMagField << std::endl;


			i=0;
			printf("Printed from host : \n");
			for (const auto& superClustersToken : superClustersTokens_) {
				for (auto& superClusRef : event.get(superClustersToken)) {
					printf("For SC i=%d Energy is :%f , theta is :%f,  \n",i,superClusRef->energy(),superClusRef->seed()->position().theta()) ;
					view[i].scSeedTheta() =  superClusRef->seed()->position().theta();
					view[i].scPhi() = superClusRef->position().phi();
					view[i].scR() = superClusRef->position().r();
					view[i].scEnergy() = superClusRef->energy();
					i++;
				}
			}	

			alpaka::memcpy(event.queue(), deviceProduct.buffer(), hostProduct.buffer());

			// Print the SoA 
			algo_.print(event.queue(), deviceProduct);
			//algo_.matchSeeds(event.queue(), deviceProduct,tracks_d.view());

			event.emplace(deviceToken_, std::move(deviceProduct));

		}

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
			edm::ParameterSetDescription desc;
			desc.add<int32_t>("size");
			desc.add<std::vector<edm::InputTag>>("getsuperclus");
			descriptions.addWithDefaultLabel(desc);
		}

	private:

		const device::EDPutToken<reco::SuperclusterDeviceCollection> deviceToken_;

		const int32_t size_;
		edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
		std::vector<edm::EDGetTokenT<std::vector<reco::SuperClusterRef>>> superClustersTokens_;
		//edm::EDGetTokenT<std::vector<reco::SuperCluster>> superClustersTokens_;
		// Try and print out the device SoA elements
		SuperclusterAlgo const algo_{};

  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(SuperClusterSoAProducer);
