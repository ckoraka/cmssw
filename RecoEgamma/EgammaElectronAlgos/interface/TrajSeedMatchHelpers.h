#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <vector>

namespace egamma {
  template<size_t nMaxEtaBins>
  struct MatchingCuts {
    using array_type = std::array<double, nMaxEtaBins>;

  public:
    explicit MatchingCuts(const array_type& dPhiHighEt,
                          const array_type& dPhiHighEtThres,
                          const array_type& dPhiLowEtGrad,
                          const array_type& dRZHighEt,
                          const array_type& dRZHighEtThres,
                          const array_type& dRZLowEtGrad,
                          const array_type& etaBins,
                          const uint nEtaBins)
      : dPhiHighEt_(dPhiHighEt),
        dPhiHighEtThres_(dPhiHighEtThres),
        dPhiLowEtGrad_(dPhiLowEtGrad),
        dRZHighEt_(dRZHighEt),
        dRZHighEtThres_(dRZHighEtThres),
        dRZLowEtGrad_(dRZLowEtGrad),
        etaBins_(etaBins),
        nEtaBins_(nEtaBins)
        {}

    virtual ~MatchingCuts() {}
    // bool operator()(const SCHitMatch& scHitMatch) const{
    //   size_t binNr = getBinNr(scHitMatch.eta);
    //   float dPhiMax = getCutValue(scHitMatch.et, dPhiHighEt_[binNr], dPhiHighEtThres_[binNr], dPhiLowEtGrad_[binNr]);
    //   if (dPhiMax >= 0 && std::abs(scHitMatch.dPhi) > dPhiMax) {
    //     return false;
    //   }
    //   float dRZMax = getCutValue(scHitMatch.et, dRZHighEt_[binNr], dRZHighEtThres_[binNr], dRZLowEtGrad_[binNr]);
    //   if (dRZMax >= 0 && std::abs(scHitMatch.dRZ) > dRZMax) {
    //     return false;
    //   }
    //   return true;
    // };

  private:
    // size_t getBinNr(float eta) const{
    //   const float absEta = std::abs(eta);
    //   for (size_t etaNr = 0; etaNr < nEtaBins_; etaNr++) {
    //     if (absEta < etaBins_[etaNr])
    //       return etaNr;
    //   }
    //   return nEtaBins_;
    // }

    // float getCutValue(float et, float highEt, float highEtThres, float lowEtGrad) const {
    //   return highEt + std::min(0.f, et - highEtThres) * lowEtGrad;
    // }

    array_type dPhiHighEt_, dPhiHighEtThres_, dPhiLowEtGrad_;
    array_type dRZHighEt_, dRZHighEtThres_, dRZLowEtGrad_;
    array_type etaBins_;
    uint nEtaBins_;
  };

  edm::ParameterSetDescription makeMatchingCugsPSetDescription() {
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
  template<size_t nMaxEtaBins>
  auto makeMatchingCuts(std::vector<edm::ParameterSet> const& cutsPSets) {
    std::vector<std::unique_ptr<egamma::MatchingCuts<nMaxEtaBins>>> matchingCuts;

    for (const auto& cutPSet : cutsPSets) {
      const auto& dPhiHighEt = cutPSet.getParameter<std::vector<double> >("dPhiMaxHighEt");
      const auto& dPhiHighEtThres = cutPSet.getParameter<std::vector<double> >("dPhiMaxHighEtThres");
      const auto& dPhiLowEtGrad = cutPSet.getParameter<std::vector<double> >("dPhiMaxLowEtGrad");
      const auto& dRZHighEt = cutPSet.getParameter<std::vector<double> >("dRZMaxHighEt");
      const auto& dRZHighEtThres = cutPSet.getParameter<std::vector<double> >("dRZMaxHighEtThres");
      const auto& dRZLowEtGrad = cutPSet.getParameter<std::vector<double> >("dRZMaxLowEtGrad");
      const auto& etaBins = cutPSet.getParameter<std::vector<double> >("etaBins");

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
      binSizeCheck(etaBins.size(), dPhiHighEt, "dPhiMaxHighEt");
      binSizeCheck(etaBins.size(), dPhiHighEtThres, "dPhiMaxHighEtThres");
      binSizeCheck(etaBins.size(), dPhiLowEtGrad, "dPhiMaxLowEtGrad");
      binSizeCheck(etaBins.size(), dRZHighEt, "dRZMaxHighEt");
      binSizeCheck(etaBins.size(), dRZHighEtThres, "dRZMaxHighEtThres");
      binSizeCheck(etaBins.size(), dRZLowEtGrad, "dRZMaxLowEtGrad");

      using array_type = std::array<double, nMaxEtaBins>;
      array_type adPhiHighEt; adPhiHighEt.fill(0);
      array_type adPhiHighEtThres; adPhiHighEtThres.fill(0);
      array_type adPhiLowEtGrad; adPhiLowEtGrad.fill(0);
      array_type adRZHighEt; adRZHighEt.fill(0);
      array_type adRZHighEtThres; adRZHighEtThres.fill(0);
      array_type adRZLowEtGrad; adRZLowEtGrad.fill(0);
      array_type aetaBins; aetaBins.fill(0);

      for (uint bin = 0; bin < etaBins.size(); bin++) {
        aetaBins[bin] = etaBins[bin];
      }

      for (uint bin = 0; bin <= etaBins.size(); bin++) {
        adPhiHighEt[bin]      = dPhiHighEt[bin];
        adPhiHighEtThres[bin] = dPhiHighEtThres[bin];
        adPhiLowEtGrad[bin]   = dPhiLowEtGrad[bin];
        adRZHighEt[bin]       = dRZHighEt[bin];
        adRZHighEtThres[bin]  = dRZHighEtThres[bin];
        adRZLowEtGrad[bin]    = dRZLowEtGrad[bin];
      }

      matchingCuts.emplace_back(std::make_unique<egamma::MatchingCuts<nMaxEtaBins>>(
        adPhiHighEt, adPhiHighEtThres, adPhiLowEtGrad,
        adRZHighEt, adRZHighEtThres, adRZLowEtGrad,
        aetaBins,
        etaBins.size()
      ));
    }

    return matchingCuts;
  }
}