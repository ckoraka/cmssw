// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/alpaka/MagneticFieldParabolicPortable.h"

//#include "RecoEgamma/EgammaElectronProducers/interface/helixToBarrelPropagator.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/alpaka/helixPropagator.h"

#include "SuperclusterAlgo.h"

using Vector3f = Eigen::Matrix<float, 3, 1>;

static constexpr float epsilon = 0.001;
/** Hard-wired numbers defining the surfaces on which the crystal front faces lie. */
static constexpr float barrelRadius = 129.f;       // p81, p50, ECAL TDR
static constexpr float barrelHalfLength = 270.9f;  // p81, p50, ECAL TDR
static constexpr float endcapRadius = 171.1f;      // fig 3.26, p81, ECAL TDR
static constexpr float endcapZ = 320.5f;           // fig 3.26, p81, ECAL TDR

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

	//--- Kernel for printing the SC SoA
  class SuperclusterAlgoKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  reco::SuperclusterDeviceCollection::View view,
                                  int32_t size) const 
		{
			const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
			// Make a strided loop over the kernel grid, covering up to "size" elements
			printf("Printed from device : \n");

			for (int32_t i : elements_with_stride(acc, size)) 
			{
				printf("For SC i=%d Energy is :%f , theta is :%f , and r is %lf and phi is %lf,  \n",i,view[i].scEnergy(),view[i].scSeedTheta(), view[i].scR(),view[i].scPhi() ) ;
				float x = view[i].scR() * alpaka::math::sin(acc,view[i].scSeedTheta()) * alpaka::math::sin(acc,view[i].scPhi());
				float y = view[i].scR() * alpaka::math::sin(acc,view[i].scSeedTheta()) * alpaka::math::sin(acc,view[i].scPhi());
				float z = view[i].scR() * alpaka::math::sin(acc,view[i].scSeedTheta());
				printf("x: %lf,  y: %lf,  z %lf",x,z,y);
				Vector3f position{x,y,z};
				printf("Value of perp2 %lf",x*x+y*y);
				printf("Calculate the magnetic field with the parabolic approx. at the SC position : %f\n",MagneticFieldParabolicPortable::MagneticFieldAtPoint(position));
			}
		}
	};


	//--- Kernel for performing the seed to SC match
	//--- Still dummy and just testing out some things
	class SeedToSuperClusterMatcher {
	public:
		template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
		ALPAKA_FN_ACC void operator()(TAcc const& acc,
										reco::SuperclusterDeviceCollection::View viewSCs,
										int32_t sizeSCs) const 
		{

			const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];

			printf("Printing from the SeedToSuperClusterMatcher kernel :");



			Vector3f x2 = {1,1,1}; // position
			Vector3f p2 = {1,1,1}; //momentun
			Vector3f gp = {0.,0.,0.,}; //position 
			Vector3f gm = {1.,1.,1.,}; //momentun 

      		double s=0; //path length
			double rho = 1; // random now for testing
			bool theSolExists = true; // random now for testing
      		Propagators::helixBarrelCylinderCrossing(gp,gm,rho,barrelRadius,theSolExists,x2,p2,s);
			printf("Position : %lf and direction : %lf and path length : %lf",x2(0),p2(0),s);


			// Strided loop over the seeds 
			//for(int32_t j : elements_with_stride(acc, viewSCs.metadata().size())) 
			//{


			//}
		}
	};

	//---- Kernel launch for printing the SC SoA collection
	void SuperclusterAlgo::print(Queue& queue, reco::SuperclusterDeviceCollection& collection) const {
		uint32_t items = 32;
		uint32_t groups = divide_up_by(collection->metadata().size(), items);
		auto workDiv = make_workdiv<Acc1D>(groups, items);
		alpaka::exec<Acc1D>(queue, workDiv, SuperclusterAlgoKernel{}, collection.view(), collection->metadata().size());
	}


	//---- Kernel launch for SC and seed matching
	void SuperclusterAlgo::matchSeeds(Queue& queue, reco::SuperclusterDeviceCollection& collection) const 
	{
		uint32_t items = 32;
		uint32_t groups = divide_up_by(collection->metadata().size(), items);
		auto workDiv = make_workdiv<Acc1D>(groups, items);
		alpaka::exec<Acc1D>(queue, workDiv, SeedToSuperClusterMatcher{},collection.view(),collection->metadata().size());
	}

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
