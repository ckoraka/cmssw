#ifndef RecoEgamma_EgammaElectronAlgos_TrajSeedMatchHelpers_h
#define RecoEgamma_EgammaElectronAlgos_TrajSeedMatchHelpers_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/TrajSeedMatchingCutsPortable.h"

#include <vector>
#include <array>

namespace egamma {

  inline edm::ParameterSetDescription makeMatchingCugsPSetDescription() {
    edm::ParameterSetDescription desc;

    edm::ParameterSetDescription cutsDesc;
    cutsDesc.add<std::vector<double> >("dPhiMaxHighEt",       {0.003});
    cutsDesc.add<std::vector<double> >("dPhiMaxHighEtThres",  {0.0});
    cutsDesc.add<std::vector<double> >("dPhiMaxLowEtGrad",    {0.0});
    cutsDesc.add<std::vector<double> >("dRZMaxHighEt",        {0.005});
    cutsDesc.add<std::vector<double> >("dRZMaxHighEtThres",   {30});
    cutsDesc.add<std::vector<double> >("dRZMaxLowEtGrad",     {-0.002});
    cutsDesc.add<std::vector<double> >("etaBins",             {});

    edm::ParameterSet defaults;
    defaults.addParameter<std::vector<double> >("dPhiMaxHighEt",      std::vector<double>{0.003});
    defaults.addParameter<std::vector<double> >("dPhiMaxHighEtThres", std::vector<double>{0.0});
    defaults.addParameter<std::vector<double> >("dPhiMaxLowEtGrad",   std::vector<double>{0.0});
    defaults.addParameter<std::vector<double> >("dRZMaxHighEt",       std::vector<double>{0.005});
    defaults.addParameter<std::vector<double> >("dRZMaxHighEtThres",  std::vector<double>{30});
    defaults.addParameter<std::vector<double> >("dRZMaxLowEtGrad",    std::vector<double>{-0.002});
    defaults.addParameter<std::vector<double> >("etaBins",            std::vector<double>{});

    desc.addVPSet("matchingCuts", cutsDesc, std::vector<edm::ParameterSet>{defaults, defaults, defaults});
    return desc;
  }

}


namespace {
  template<size_t nMaxEtaBins, size_t nMaxHits>
  auto makeMatchingCuts(std::vector<edm::ParameterSet> const& cutsPSets) {
    std::array<egamma::MatchingCuts<nMaxEtaBins>, nMaxHits> matchingCuts;

    if (cutsPSets.size() != nMaxHits) {
      throw cms::Exception("InvalidConfig")
          << " when constructing egamma::MatchingCuts " << cutsPSets.size() << " ParameterSets were provided, "
          << "but only nMaxHits " << nMaxHits << "are expected.";
    }

    for (size_t nHit = 0; nHit < nMaxHits; nHit++) {
      const auto& cutPSet = cutsPSets[nHit];
      const auto& vdPhiHighEt = cutPSet.getParameter<std::vector<double> >("dPhiMaxHighEt");
      const auto& vdPhiHighEtThres = cutPSet.getParameter<std::vector<double> >("dPhiMaxHighEtThres");
      const auto& vdPhiLowEtGrad = cutPSet.getParameter<std::vector<double> >("dPhiMaxLowEtGrad");
      const auto& vdRZHighEt = cutPSet.getParameter<std::vector<double> >("dRZMaxHighEt");
      const auto& vdRZHighEtThres = cutPSet.getParameter<std::vector<double> >("dRZMaxHighEtThres");
      const auto& vdRZLowEtGrad = cutPSet.getParameter<std::vector<double> >("dRZMaxLowEtGrad");
      const auto& vetaBins = cutPSet.getParameter<std::vector<double> >("etaBins");

      auto binSizeCheck = [](size_t sizeEtaBins, const std::vector<double>& vec, const std::string& name) {
        if (vec.size() != sizeEtaBins + 1) {
          throw cms::Exception("InvalidConfig")
              << " when constructing egamma::MatchingCuts " << name << " has " << vec.size()
              << " bins, it should be equal to #bins of etaBins+1 " << sizeEtaBins + 1;
        }
        if (sizeEtaBins + 1 > nMaxEtaBins) {
          throw cms::Exception("InvalidConfig")
              << " when constructing egamma::MatchingCuts " << name << " has " << sizeEtaBins + 1
              << " bins, it should be at most nMaxEtaBins " << nMaxEtaBins;

        }
      };
      const auto nvetaBins = vetaBins.size();
      binSizeCheck(nvetaBins, vdPhiHighEt, "dPhiMaxHighEt");
      binSizeCheck(nvetaBins, vdPhiHighEtThres, "dPhiMaxHighEtThres");
      binSizeCheck(nvetaBins, vdPhiLowEtGrad, "dPhiMaxLowEtGrad");
      binSizeCheck(nvetaBins, vdRZHighEt, "dRZMaxHighEt");
      binSizeCheck(nvetaBins, vdRZHighEtThres, "dRZMaxHighEtThres");
      binSizeCheck(nvetaBins, vdRZLowEtGrad, "dRZMaxLowEtGrad");

      egamma::MatchingCuts<nMaxEtaBins> cuts{};
      cuts.nEtaBins_ = nvetaBins;

      for (uint bin = 0; bin < nvetaBins; bin++) {
        cuts.etaBins_[bin] = vetaBins[bin];
      }

      for (uint bin = 0; bin <= nvetaBins; bin++) {
        cuts.dPhiHighEt_[bin]      = vdPhiHighEt[bin];
        cuts.dPhiHighEtThres_[bin] = vdPhiHighEtThres[bin];
        cuts.dPhiLowEtGrad_[bin]   = vdPhiLowEtGrad[bin];
        cuts.dRZHighEt_[bin]       = vdRZHighEt[bin];
        cuts.dRZHighEtThres_[bin]  = vdRZHighEtThres[bin];
        cuts.dRZLowEtGrad_[bin]    = vdRZLowEtGrad[bin];
      }

      matchingCuts[nHit] = cuts;
    }

    return matchingCuts;
  }
}

#endif
