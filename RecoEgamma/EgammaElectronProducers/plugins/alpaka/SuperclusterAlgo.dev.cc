// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <alpaka/alpaka.hpp>

#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "SuperclusterAlgo.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

  class SuperclusterAlgoKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  portableSuperclusterSoA::SuperclusterDeviceCollection::View view,
                                  int32_t size) const {

      const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
      // make a strided loop over the kernel grid, covering up to "size" elements
      printf("Printed from device : ");
      for (int32_t i : elements_with_stride(acc, size)) {
        printf("For SC i=%d Energy is :%f , theta is :%f,  \n",i,view[i].scEnergy(),view[i].scSeedTheta()) ;
      }
    }
  };

  void SuperclusterAlgo::print(Queue& queue, portableSuperclusterSoA::SuperclusterDeviceCollection& collection) const {
    uint32_t items = 32;
    uint32_t groups = divide_up_by(collection->metadata().size(), items);
    auto workDiv = make_workdiv<Acc1D>(groups, items);
    alpaka::exec<Acc1D>(queue, workDiv, SuperclusterAlgoKernel{}, collection.view(), collection->metadata().size());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
