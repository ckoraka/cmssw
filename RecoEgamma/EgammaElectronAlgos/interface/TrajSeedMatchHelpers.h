#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <vector>

namespace egamma {
  class MatchingCuts {
  public:
    explicit MatchingCuts(const edm::ParameterSet& pset) 
      : dPhiHighEt_(pset.getParameter<std::vector<double> >("dPhiMaxHighEt")),
        dPhiHighEtThres_(pset.getParameter<std::vector<double> >("dPhiMaxHighEtThres")),
        dPhiLowEtGrad_(pset.getParameter<std::vector<double> >("dPhiMaxLowEtGrad")),
        dRZHighEt_(pset.getParameter<std::vector<double> >("dRZMaxHighEt")),
        dRZHighEtThres_(pset.getParameter<std::vector<double> >("dRZMaxHighEtThres")),
        dRZLowEtGrad_(pset.getParameter<std::vector<double> >("dRZMaxLowEtGrad")),
        etaBins_(pset.getParameter<std::vector<double> >("etaBins")) {}

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
    //   for (size_t etaNr = 0; etaNr < etaBins_.size(); etaNr++) {
    //     if (absEta < etaBins_[etaNr])
    //       return etaNr;
    //   }
    //   return etaBins_.size();
    // }

    // float getCutValue(float et, float highEt, float highEtThres, float lowEtGrad) const {
    //   return highEt + std::min(0.f, et - highEtThres) * lowEtGrad;
    // }

  private:
    std::vector<double> dPhiHighEt_, dPhiHighEtThres_, dPhiLowEtGrad_;
    std::vector<double> dRZHighEt_, dRZHighEtThres_, dRZLowEtGrad_;
    std::vector<double> etaBins_;
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
  auto makeMatchingCuts(std::vector<edm::ParameterSet> const& cutsPSets) {
    std::vector<std::unique_ptr<egamma::MatchingCuts> > matchingCuts;

    for (const auto& cutPSet : cutsPSets) {
      matchingCuts.emplace_back(std::make_unique<egamma::MatchingCuts>(cutPSet));
    }

    return matchingCuts;
  }
}