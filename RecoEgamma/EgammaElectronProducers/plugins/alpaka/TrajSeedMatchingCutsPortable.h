#ifndef RecoEgamma_EgammaElectronProducers_plugins_alpaka_TrajSeedMatchingCutsPortable_h
#define RecoEgamma_EgammaElectronProducers_plugins_alpaka_TrajSeedMatchingCutsPortable_h

#include <array>

namespace egamma {
  template<size_t nMaxEtaBins>
  struct MatchingCuts {
    using array_type = std::array<double, nMaxEtaBins>;

  public:
    MatchingCuts() {};

    array_type dPhiHighEt_{}, dPhiHighEtThres_{}, dPhiLowEtGrad_{};
    array_type dRZHighEt_{}, dRZHighEtThres_{}, dRZLowEtGrad_{};
    array_type etaBins_{};
    uint nEtaBins_{};
  };

}

#endif
