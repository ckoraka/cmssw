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

  private:
    // size_t getBinNr(float eta) const{
    //   const float absEta = std::abs(eta);
    //   for (size_t etaNr = 0; etaNr < nEtaBins_; etaNr++) {
    //     if (absEta < etaBins_[etaNr])
    //       return etaNr;
    //   }
    //   return nEtaBins_;
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

      using array_type = std::array<double, nMaxEtaBins>;
      array_type dPhiHighEt;      dPhiHighEt.fill(0);
      array_type dPhiHighEtThres; dPhiHighEtThres.fill(0);
      array_type dPhiLowEtGrad;   dPhiLowEtGrad.fill(0);
      array_type dRZHighEt;       dRZHighEt.fill(0);
      array_type dRZHighEtThres;  dRZHighEtThres.fill(0);
      array_type dRZLowEtGrad;    dRZLowEtGrad.fill(0);
      array_type etaBins;         etaBins.fill(0);

      for (uint bin = 0; bin < nvetaBins; bin++) {
        etaBins[bin] = vetaBins[bin];
      }

      for (uint bin = 0; bin <= nvetaBins; bin++) {
        dPhiHighEt[bin]      = vdPhiHighEt[bin];
        dPhiHighEtThres[bin] = vdPhiHighEtThres[bin];
        dPhiLowEtGrad[bin]   = vdPhiLowEtGrad[bin];
        dRZHighEt[bin]       = vdRZHighEt[bin];
        dRZHighEtThres[bin]  = vdRZHighEtThres[bin];
        dRZLowEtGrad[bin]    = vdRZLowEtGrad[bin];
      }

      matchingCuts.emplace_back(std::make_unique<egamma::MatchingCuts<nMaxEtaBins>>(
        dPhiHighEt, dPhiHighEtThres, dPhiLowEtGrad,
        dRZHighEt, dRZHighEtThres, dRZLowEtGrad,
        etaBins,
        nvetaBins
      ));
    }

    return matchingCuts;
  }
}