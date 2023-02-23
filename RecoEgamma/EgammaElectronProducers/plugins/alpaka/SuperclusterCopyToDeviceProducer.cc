#include "DataFormats/Portable/interface/Product.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterSoA.h"
#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/SuperclusterHostCollection.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


#include "HeterogeneousCore/AlpakaCore/interface/ScopedContext.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaServices/interface/alpaka/AlpakaService.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/global/EDProducer.h"

#include "SuperclusterAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

	class SuperclusterCopyToDeviceProducer : public global::EDProducer<> {
	public:
		SuperclusterCopyToDeviceProducer(edm::ParameterSet const& config)
			: deviceToken_{produces()},
			size_{config.getParameter<int32_t>("size")}{
			//superClustersTokens_ = consumes(config.getParameter<edm::InputTag>("getsuperclus"));
			for (const auto& scTag : config.getParameter<std::vector<edm::InputTag>>("getsuperclus")) {
				superClustersTokens_.emplace_back(consumes(scTag));
			}
		}
  
		void produce(edm::StreamID sid, device::Event& event, device::EventSetup const&) const override {

			int i=0;
			printf("Printed from host : \n");
			for (const auto& superClustersToken : superClustersTokens_) {
				for (auto& superClusRef : event.get(superClustersToken)) {
					i++;
				}
			}	
			
			portableSuperclusterSoA::SuperclusterHostCollection hostProduct{i, event.queue()};
			portableSuperclusterSoA::SuperclusterDeviceCollection deviceProduct{i, event.queue()};

			auto& view = hostProduct.view();

			/*
			int i=0;
			for (auto& superClusRef : event.get(superClustersTokens_)) {
				view[i].scSeedTheta() =  superClusRef.seed()->position().theta();
				view[i].scPhi() = superClusRef.position().phi();
				view[i].scR() = superClusRef.position().r();
				view[i].scEnergy() = superClusRef.energy();
				i++;
			}
			*/

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

			// run the algorithm
			algo_.print(event.queue(), deviceProduct);

			event.emplace(deviceToken_, std::move(deviceProduct));

		}

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
			edm::ParameterSetDescription desc;
			desc.add<int32_t>("size");
			desc.add<std::vector<edm::InputTag>>("getsuperclus");
			descriptions.addWithDefaultLabel(desc);
		}

	private:

		const device::EDPutToken<portableSuperclusterSoA::SuperclusterDeviceCollection> deviceToken_;

		const int32_t size_;
		std::vector<edm::EDGetTokenT<std::vector<reco::SuperClusterRef>>> superClustersTokens_;
		//edm::EDGetTokenT<std::vector<reco::SuperCluster>> superClustersTokens_;

		// Try and print out the device SoA elements
		SuperclusterAlgo const algo_{};
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#include "HeterogeneousCore/AlpakaCore/interface/MakerMacros.h"
DEFINE_FWK_ALPAKA_MODULE(SuperclusterCopyToDeviceProducer);
