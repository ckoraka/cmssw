#ifndef DataFormats_Track_PixelTrackUtilities_h
#define DataFormats_Track_PixelTrackUtilities_h

#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "DataFormats/Track/interface/PixelTrackDefinitions.h"
#include "DataFormats/Track/interface/PixelTrackLayout.h"

//Copied from Breno's branch

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  // Methods that operate on View and ConstView of the TrackSoA, and cannot be class methods.
  template <typename TrackerTraits>
  struct TracksUtilities {
    using TrackSoAView = typename TrackSoA<TrackerTraits>::template TrackSoAHeterogeneousLayout<>::View;
    using TrackSoAConstView = typename TrackSoA<TrackerTraits>::template TrackSoAHeterogeneousLayout<>::ConstView;
    using hindex_type = typename TrackSoA<TrackerTraits>::hindex_type;

    // State at the Beam spot
    // phi,tip,1/pt,cotan(theta),zip
    static ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE float charge(const TrackSoAConstView &tracks, int32_t i) {
      return std::copysign(1.f, tracks[i].state()(2));
    }

    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE float phi(const TrackSoAConstView &tracks, int32_t i) {
      return tracks[i].state()(0);
    }

    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE float tip(const TrackSoAConstView &tracks, int32_t i) {
      return tracks[i].state()(1);
    }

    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE float zip(const TrackSoAConstView &tracks, int32_t i) {
      return tracks[i].state()(4);
    }

    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE bool isTriplet(const TrackSoAConstView &tracks, int i) {
      return tracks[i].nLayers() == 3;
    }

    template <typename V3, typename M3, typename V2, typename M2>
    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE void copyFromCircle(
        TrackSoAView &tracks, V3 const &cp, M3 const &ccov, V2 const &lp, M2 const &lcov, float b, int32_t i) {
      tracks[i].state() << cp.template cast<float>(), lp.template cast<float>();

      tracks[i].state()(2) = tracks[i].state()(2) * b;
      auto cov = tracks[i].covariance();
      cov(0) = ccov(0, 0);
      cov(1) = ccov(0, 1);
      cov(2) = b * float(ccov(0, 2));
      cov(4) = cov(3) = 0;
      cov(5) = ccov(1, 1);
      cov(6) = b * float(ccov(1, 2));
      cov(8) = cov(7) = 0;
      cov(9) = b * b * float(ccov(2, 2));
      cov(11) = cov(10) = 0;
      cov(12) = lcov(0, 0);
      cov(13) = lcov(0, 1);
      cov(14) = lcov(1, 1);
    }

    template <typename V5, typename M5>
    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE void copyFromDense(TrackSoAView &tracks,
                                                                            V5 const &v,
                                                                            M5 const &cov,
                                                                            int32_t i) {
      tracks[i].state() = v.template cast<float>();
      for (int j = 0, ind = 0; j < 5; ++j)
        for (auto k = j; k < 5; ++k)
          tracks[i].covariance()(ind++) = cov(j, k);
    }

    template <typename V5, typename M5>
    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE void copyToDense(const TrackSoAConstView &tracks,
                                                                          V5 &v,
                                                                          M5 &cov,
                                                                          int32_t i) {
      v = tracks[i].state().template cast<typename V5::Scalar>();
      for (int j = 0, ind = 0; j < 5; ++j) {
        cov(j, j) = tracks[i].covariance()(ind++);
        for (auto k = j + 1; k < 5; ++k)
          cov(k, j) = cov(j, k) = tracks[i].covariance()(ind++);
      }
    }

    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE int computeNumberOfLayers(const TrackSoAConstView &tracks,
                                                                                   int32_t i) {
      auto pdet = tracks.detIndices().begin(i);
      int nl = 1;
      auto ol = pixelTopology::getLayer<TrackerTraits>(*pdet);
      for (; pdet < tracks.detIndices().end(i); ++pdet) {
        auto il = pixelTopology::getLayer<TrackerTraits>(*pdet);
        if (il != ol)
          ++nl;
        ol = il;
      }
      return nl;
    }

    static constexpr ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE int nHits(const TrackSoAConstView &tracks, int i) {
      return tracks.detIndices().size(i);
    }
  };

  namespace pixelTrack {

    template <typename TrackerTraits, typename Enable = void>
    struct QualityCutsT {};

    template <typename TrackerTraits>
    struct QualityCutsT<TrackerTraits, pixelTopology::isPhase1Topology<TrackerTraits>> {
      using TrackSoAView = typename TrackSoA<TrackerTraits>::template TrackSoAHeterogeneousLayout<>::View;
      using TrackSoAConstView = typename TrackSoA<TrackerTraits>::template TrackSoAHeterogeneousLayout<>::ConstView;
      using tracksHelper = TracksUtilities<TrackerTraits>;
      // chi2 cut = chi2Scale * (chi2Coeff[0] + pT/GeV * (chi2Coeff[1] + pT/GeV * (chi2Coeff[2] + pT/GeV * chi2Coeff[3])))
      float chi2Coeff[4];
      float chi2MaxPt;  // GeV
      float chi2Scale;

      struct Region {
        float maxTip;  // cm
        float minPt;   // GeV
        float maxZip;  // cm
      };

      Region triplet;
      Region quadruplet;

      ALPAKA_FN_ACC ALPAKA_FN_INLINE bool isHP(const TrackSoAConstView &tracks, int nHits, int it) const {
        // impose "region cuts" based on the fit results (phi, Tip, pt, cotan(theta)), Zip)
        // default cuts:
        //   - for triplets:    |Tip| < 0.3 cm, pT > 0.5 GeV, |Zip| < 12.0 cm
        //   - for quadruplets: |Tip| < 0.5 cm, pT > 0.3 GeV, |Zip| < 12.0 cm
        // (see CAHitNtupletGeneratorGPU.cc)
        auto const &region = (nHits > 3) ? quadruplet : triplet;
        return (std::abs(tracksHelper::tip(tracks, it)) < region.maxTip) and (tracks.pt(it) > region.minPt) and
               (std::abs(tracksHelper::zip(tracks, it)) < region.maxZip);
      }

      ALPAKA_FN_ACC ALPAKA_FN_INLINE bool strictCut(const TrackSoAConstView &tracks, int it) const {
        auto roughLog = [](float x) {
          // max diff [0.5,12] at 1.25 0.16143
          // average diff  0.0662998
          union IF {
            uint32_t i;
            float f;
          };
          IF z;
          z.f = x;
          uint32_t lsb = 1 < 21;
          z.i += lsb;
          z.i >>= 21;
          auto f = z.i & 3;
          int ex = int(z.i >> 2) - 127;

          // log2(1+0.25*f)
          // averaged over bins
          const float frac[4] = {0.160497f, 0.452172f, 0.694562f, 0.901964f};
          return float(ex) + frac[f];
        };

        float pt = std::min<float>(tracks.pt(it), chi2MaxPt);
        float chi2Cut = chi2Scale * (chi2Coeff[0] + roughLog(pt) * chi2Coeff[1]);
        if (tracks.chi2(it) >= chi2Cut) {
#ifdef NTUPLE_FIT_DEBUG
          printf("Bad chi2 %d pt %f eta %f chi2 %f\n", it, tracks.pt(it), tracks.eta(it), tracks.chi2(it));
#endif
          return true;
        }
        return false;
      }
    };

    template <typename TrackerTraits>
    struct QualityCutsT<TrackerTraits, pixelTopology::isPhase2Topology<TrackerTraits>> {
      using TrackSoAView = typename TrackSoA<TrackerTraits>::template TrackSoAHeterogeneousLayout<>::View;
      using TrackSoAConstView = typename TrackSoA<TrackerTraits>::template TrackSoAHeterogeneousLayout<>::ConstView;
      using tracksHelper = TracksUtilities<TrackerTraits>;

      float maxChi2;
      float minPt;
      float maxTip;
      float maxZip;

      ALPAKA_FN_ACC ALPAKA_FN_INLINE bool isHP(const TrackSoAConstView &tracks, int nHits, int it) const {
        return (std::abs(tracksHelper::tip(tracks, it)) < maxTip) and (tracks.pt(it) > minPt) and
               (std::abs(tracksHelper::zip(tracks, it)) < maxZip);
      }
      ALPAKA_FN_ACC ALPAKA_FN_INLINE bool strictCut(const TrackSoAConstView &tracks, int it) const {
        return tracks.chi2(it) >= maxChi2;
      }
    };

  }  // namespace pixelTrack

  // TODO: Should those be placed in the ALPAKA_ACCELERATOR_NAMESPACE
  template struct TracksUtilities<pixelTopology::Phase1>;
  template struct TracksUtilities<pixelTopology::Phase2>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
