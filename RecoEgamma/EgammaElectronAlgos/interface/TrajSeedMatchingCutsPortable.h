#ifndef RecoEgamma_EgammaElectronAlgos_TrajSeedMatchingCutsPortable_h
#define RecoEgamma_EgammaElectronAlgos_TrajSeedMatchingCutsPortable_h

#include <array>

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

    array_type dPhiHighEt_{}, dPhiHighEtThres_{}, dPhiLowEtGrad_{};
    array_type dRZHighEt_{}, dRZHighEtThres_{}, dRZLowEtGrad_{};
    array_type etaBins_{};
    uint nEtaBins_{};
  };

}

#endif
