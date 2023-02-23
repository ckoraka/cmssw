#ifndef RecoEgamma_EgammaElectronProducers_plugins_alpaka_SuperclusterAlgo_h
#define RecoEgamma_EgammaElectronProducers_plugins_alpaka_SuperclusterAlgo_h

#include <Eigen/Core>

#include "DataFormats/PortableTestObjects/interface/alpaka/TestDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "DataFormats/Track/interface/PixelTrackDefinitions.h"
#include "DataFormats/Track/interface/TrackSoAHost.h" 
#include "DataFormats/Track/interface/alpaka/TrackSoADevice.h" 


namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class SuperclusterAlgo {
  public:
    void print(Queue& queue, portableSuperclusterSoA::SuperclusterDeviceCollection& collection) const;
  };

	class FillTrackSoA {
	public:
	  void fillTrackSoA(TrackSoAView<pixelTopology::Phase1> tracks_view, Queue& queue) ;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif 
