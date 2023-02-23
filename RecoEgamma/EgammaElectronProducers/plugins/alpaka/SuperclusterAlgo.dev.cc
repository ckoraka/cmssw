// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

// Includes to create and fill the dummy track SoA
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "DataFormats/Track/interface/PixelTrackDefinitions.h"
#include "DataFormats/Track/interface/alpaka/TrackSoADevice.h"
#include "DataFormats/Track/interface/TrackSoAHost.h" 

#include "SuperclusterAlgo.h"

using Quality = pixelTrack::Quality;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

	//--- Kernel for printing the SC SoA
  class SuperclusterAlgoKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  portableSuperclusterSoA::SuperclusterDeviceCollection::View view,
                                  int32_t size) const {

      const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
      // make a strided loop over the kernel grid, covering up to "size" elements
      printf("Printed from device : \n");
      for (int32_t i : elements_with_stride(acc, size)) {
        printf("For SC i=%d Energy is :%f , theta is :%f,  \n",i,view[i].scEnergy(),view[i].scSeedTheta()) ;
      }
    }
  };

	//--- Kernel for performing the seed to SC match
	class SeedToSuperClusterMatcher {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  TrackSoAView<pixelTopology::Phase1> view,
                                  int32_t size,
                                  portableSuperclusterSoA::SuperclusterDeviceCollection::View viewSCs,
                                  int32_t sizeSCs) const {

      const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];

			printf("Printing from the SeedToSuperClusterMatcher kernel :");
			// Strided loop over the seeds 
			for(int32_t j : elements_with_stride(acc, view.metadata().size())) {
        printf("Track pT : %f \n",view[j].pt());
				for(int32_t i=0;i<=sizeSCs;++i){
					printf("SC Energy : %f \n",viewSCs[i].scEnergy());
					// Algo implementation should go here :
					//
				}
      }

    }
  };

	//---- Kernel launch for printing the SC SoA collection
  void SuperclusterAlgo::print(Queue& queue, portableSuperclusterSoA::SuperclusterDeviceCollection& collection) const {
    uint32_t items = 32;
    uint32_t groups = divide_up_by(collection->metadata().size(), items);
    auto workDiv = make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, SuperclusterAlgoKernel{}, collection.view(), collection->metadata().size());
  }

	//---- Kernel launch for SC and seed matching
  void SuperclusterAlgo::matchSeeds(Queue& queue, portableSuperclusterSoA::SuperclusterDeviceCollection& collection,
																		TrackSoAView<pixelTopology::Phase1> tracks_view) const {
    uint32_t items = 32;
    uint32_t groups = divide_up_by(tracks_view.metadata().size(), items);
    auto workDiv = make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, SeedToSuperClusterMatcher{}, tracks_view, tracks_view.metadata().size(),collection.view(),collection->metadata().size());
  }


	//----------- Kernel and kernel launch for track SoA
	template <typename TrackerTraits>
	class FillTrackSoAKernel {
	public:
		template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
		ALPAKA_FN_ACC void operator()(TAcc const& acc, TrackSoAView<TrackerTraits> tracks_view) const {
			const int32_t i = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];

			if (i == 0) {
				tracks_view.nTracks() = 420;
			}

			for (int32_t j : elements_with_stride(acc, tracks_view.metadata().size())) {
				tracks_view[j].pt() = (float)j;
				tracks_view[j].eta() = (float)j;
				tracks_view[j].chi2() = (float)j;
				tracks_view[j].quality() = (Quality)(j % 256);
				tracks_view[j].nLayers() = j % 128;
				tracks_view.hitIndices().off[j] = j;
			}

			//alpaka::syncBlockThreads(acc);
			//for (int32_t j : elements_with_stride(acc, tracks_view.metadata().size())) {
      //  printf("Track pT : %f \n",tracks_view[j].pt());
      //}
		}
	};

	//---- Fill the track SoA
	void FillTrackSoA::fillTrackSoA(TrackSoAView<pixelTopology::Phase1> tracks_view, Queue& queue) {
		uint32_t items = 64;
		uint32_t groups = divide_up_by(tracks_view.metadata().size(), items);
		auto workDiv = make_workdiv<Acc1D>(groups, items);
		alpaka::exec<Acc1D>(queue, workDiv, FillTrackSoAKernel<pixelTopology::Phase1>{}, tracks_view);
	}

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
