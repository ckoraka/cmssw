#include <string>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/EgammaReco/interface/EleSeedHostCollection.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"


class ElectronSeedConverter : public edm::global::EDProducer<> {
public:
  explicit ElectronSeedConverter(const edm::ParameterSet&);

  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const final;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    const edm::EDGetTokenT<TrajectorySeedCollection> initialSeedsToken_;
    const edm::EDGetTokenT<reco::EleSeedHostCollection> matchedEleSeedSoAToken_; 
	const edm::EDGetTokenT<std::vector<reco::SuperClusterRef>> superClustersTokens_;
    const edm::EDPutTokenT<reco::ElectronSeedCollection> putToken_;
};


void ElectronSeedConverter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("initialSeeds", {"hltElePixelSeedsCombined"});
    desc.add<edm::InputTag>("eleSeedsSoA", edm::InputTag{"ElectronNHitSeedAlpakaProducer"});
  	desc.add<edm::InputTag>("superClusters", {"hltEgammaSuperClustersToPixelMatch"});
    descriptions.add("electronNHitSeedConverter", desc);
}

ElectronSeedConverter::ElectronSeedConverter(const edm::ParameterSet& pset)
    :initialSeedsToken_(consumes(pset.getParameter<edm::InputTag>("initialSeeds"))),
    matchedEleSeedSoAToken_(consumes(pset.getParameter<edm::InputTag>("eleSeedsSoA"))),    
    superClustersTokens_(consumes(pset.getParameter<edm::InputTag>("superClusters"))),
    putToken_(produces()) {}

void ElectronSeedConverter::produce(edm::StreamID, edm::Event& event, const edm::EventSetup& iSetup) const {

    auto const& view = event.get(matchedEleSeedSoAToken_).const_view();

    reco::ElectronSeedCollection eleSeeds{};
    // eleSeeds->reserve(eleSeedSoAView.metadata().size()); // This is too big 

	auto const& superClusterRefs = event.get(superClustersTokens_);
	auto const& initialSeeds = event.get(initialSeedsToken_);

    for (int i = 0; i < view.metadata().size(); ++i) 
	{
        if (view[i].isMatched() > 0) {
            int matchedScID = view[i].matchedScID();
            int seedID = view[i].id();
            
			std::cout << "  matchedScID: " << view[i].matchedScID() << std::endl;
            
			if (matchedScID >= 0 && static_cast<unsigned int>(matchedScID) < superClusterRefs.size() && seedID >= 0 && static_cast<unsigned int>(seedID) < initialSeeds.size()) 
			{
            	const reco::SuperClusterRef& superClusRef = superClusterRefs[matchedScID];
				const TrajectorySeed& matchedSeed = initialSeeds[seedID];
				reco::ElectronSeed eleSeed(matchedSeed);
				reco::ElectronSeed::CaloClusterRef caloClusRef(superClusRef);
				eleSeed.setCaloCluster(caloClusRef);
				eleSeeds.emplace_back(eleSeed);

            } else {
                edm::LogWarning("ElectronSeedConverter") << "Index out of bounds for SC ID: " << matchedScID << " or Seed ID: " << seedID;
            }
    	}
	}
    std::cout << "New eleSeeds size " << eleSeeds.size() << std::endl;
    event.emplace(putToken_, std::move(eleSeeds));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronSeedConverter);
