// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/alpaka/MagneticFieldParabolicPortable.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "SuperclusterAlgo.h"

using Vector3f = Eigen::Matrix<float, 3, 1>;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

	//--- Kernel for printing the SC SoA
  class SuperclusterAlgoKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  reco::SuperclusterDeviceCollection::View view,
                                  int32_t size) const {
      const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
      // make a strided loop over the kernel grid, covering up to "size" elements
      //printf("Printed from device : \n");
      for (int32_t i : elements_with_stride(acc, size)) {
        printf("For SC i=%d Energy is :%f , theta is :%f,  \n",i,view[i].scEnergy(),view[i].scSeedTheta()) ;
		// Try and see if MagField works :
		Vector3f position{1,1,1};
		printf("Calculate the Mag Field at the SC position : %f\n",MagneticFieldParabolicPortable::MagneticFieldAtPoint(position));

      }
    }
  };

/*
	//--- Kernel for performing the seed to SC match
	class SeedToSuperClusterMatcher {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  TrackSoAView<pixelTopology::Phase1> view,
                                  int32_t size,
                                  reco::SuperclusterDeviceCollection::View viewSCs,
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
*/
	//---- Kernel launch for printing the SC SoA collection
  void SuperclusterAlgo::print(Queue& queue, reco::SuperclusterDeviceCollection& collection) const {
    uint32_t items = 32;
    uint32_t groups = divide_up_by(collection->metadata().size(), items);
    auto workDiv = make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, SuperclusterAlgoKernel{}, collection.view(), collection->metadata().size());
  }

/*
	//---- Kernel launch for SC and seed matching
  void SuperclusterAlgo::matchSeeds(Queue& queue, reco::SuperclusterDeviceCollection& collection,
																		TrackSoAView<pixelTopology::Phase1> tracks_view) const {
    uint32_t items = 32;
    uint32_t groups = divide_up_by(tracks_view.metadata().size(), items);
    auto workDiv = make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, SeedToSuperClusterMatcher{}, tracks_view, tracks_view.metadata().size(),collection.view(),collection->metadata().size());
  }
*/

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
