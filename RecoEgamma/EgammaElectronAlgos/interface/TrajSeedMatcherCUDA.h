#ifndef TrajSeedMatcherCUDA_h
#define TrajSeedMatcherCUDA_h

#include <cstddef>
#include <cstdint>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "DataFormats/DetId/interface/DetId.h"

namespace TrajSeedMatcherCUDA {
    namespace common {
        struct GlobalPoint {
            float eta;    
            float phi;          
            float theta;         
            float x;         
            float y;         
            float z;         
        };
        enum PropagationDirection {oppositeToMomentum, alongMomentum, anyDirection, invalidDirection };
    }

  struct SCHitMatch {
    const DetId detId = 0;
    const common::GlobalPoint hitPos;
    const float dRZ = std::numeric_limits<float>::max();
    const float dPhi = std::numeric_limits<float>::max();
    //const TrackingRecHit& hit;
    const float et = 0.f;
    const float eta = 0.f;
    const float phi = 0.f;
    const int charge = 0;
    const int nrClus = 0;
  };

  struct MatchInfo {
    const DetId detId;
    const float dRZPos;
    const float dRZNeg;
    const float dPhiPos;
    const float dPhiNeg;
  };

  struct SeedWithInfo {
    //const TrajectorySeed& seed;
    const std::vector<MatchInfo> matchInfos;
    const int nrValidLayers;
  };

  struct TrajStateOnSurface{


  };

}

#endif

/*
__global__ void matcher(unsigned int nSCs, EcalSC::EcalSCSoA* SCs) {
   
   // First make the trajectories for each SCs this can go outside of this kernel and be passed?
   makeTrajStateOnSurface(candPos, energy, -1);
   
 
    // then start processing the seeds 
    processSeed(seed, candPos, energy, scTrajStateOnSurfNeg);


  //next try passing these variables in once...
  const float candEta = candPos.eta();
  const float candEt = energy * std::sin(candPos.theta());
  const int charge = initialTrajState.charge();

  std::vector<SCHitMatch> matches;
  FreeTrajectoryState firstMatchFreeTraj;
  GlobalPoint prevHitPos;  // Maybe  make a structure since it is used a lot!
  GlobalPoint vertex; 
  const auto nCuts = cfg_.matchingCuts.size();
  for (size_t iHit = 0;
       matches.size() < nCuts && iHit < seed.nHits() && (cfg_.enableHitSkipping || iHit == matches.size());
       iHit++) {
    auto const& recHit = *(seed.recHits().begin() + iHit);

    if (!recHit.isValid()) {
      continue;
    }

    const bool doFirstMatch = matches.empty();

    auto const& trajState = doFirstMatch
                                ? getTrajStateFromVtx(recHit, initialTrajState, backwardPropagator_)
                                : getTrajStateFromPoint(recHit, firstMatchFreeTraj, prevHitPos, forwardPropagator_);
    if (!trajState.isValid()) {
      continue;
    }

    auto const& vtxForMatchObject = doFirstMatch ? vprim_ : vertex;
    auto match = makeSCHitMatch(vtxForMatchObject, trajState, recHit, candEt, candEta, candPos.phi(), charge, 1);

   
   
   }
*/