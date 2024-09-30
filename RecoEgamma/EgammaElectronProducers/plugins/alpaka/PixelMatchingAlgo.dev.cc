// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include <alpaka/alpaka.hpp>
#include <Eigen/Core>

#include "DataFormats/EgammaReco/interface/alpaka/EleSeedDeviceCollection.h"
#include "DataFormats/EgammaReco/interface/alpaka/SuperclusterDeviceCollection.h"
#include "MagneticField/ParametrizedEngine/interface/ParabolicParametrizedMagneticField.h"

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/traits.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

//#include "RecoEgamma/EgammaElectronAlgos/interface/alpaka/helixPropagator.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/alpaka/helixBarrelPlaneCrossingByCircle.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/alpaka/ftsFromVertexToPointPortable.h"

#include "PixelMatchingAlgo.h"


using Vector3f = Eigen::Matrix<double, 3, 1>;

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  using namespace cms::alpakatools;

   //--- Kernel for printing the SC SoA
  class printSCSoA {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  reco::SuperclusterDeviceCollection::View view,
                                  int32_t size) const 
		{
			const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
			printf("Printed from device : \n");
			for (int32_t i : uniform_elements(acc, size)) 
			{
				printf("For SC i=%d Energy is :%f , theta is :%f , and r is %lf and phi is %lf,  \n",i,view[i].scEnergy(),view[i].scSeedTheta(), view[i].scR(),view[i].scPhi() ) ;
				float x = view[i].scR() * alpaka::math::sin(acc,view[i].scSeedTheta()) * alpaka::math::sin(acc,view[i].scPhi());
				float y = view[i].scR() * alpaka::math::sin(acc,view[i].scSeedTheta()) * alpaka::math::cos(acc,view[i].scPhi());
				float z = view[i].scR() * alpaka::math::sin(acc,view[i].scSeedTheta());
				printf("x: %lf,  y: %lf,  z %lf",x,z,y);
				Vector3f position{x,y,z};
				printf("Value of perp2 %lf \n",x*x+y*y);
				printf("Calculate the magnetic field with the parabolic approximation at the SC position : %f\n",magneticFieldParabolicPortable::magneticFieldAtPoint(position));
			}
		}
	};

   //--- Kernel for printing the electron seeds SoA
	class printElectronSeedSoA {
	public:
	template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
	ALPAKA_FN_ACC void operator()(TAcc const& acc,
									reco::EleSeedDeviceCollection::View view,
									int32_t size) const 
		{
			const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];
			printf("Printed from device : \n");
			for (int32_t i : uniform_elements(acc, size)) 
			{
				printf("For electron seed %d number of hits is %d \n ",i,view[i].nHits());
				printf("Hit 0 position is (%lf , %lf , %lf ) \n", view[i].hitPosX().x(),view[i].hitPosY().x(),view[i].hitPosZ().x());
				printf("Hit 1 position is (%lf , %lf , %lf ) \n", view[i].hitPosX().y(),view[i].hitPosY().y(),view[i].hitPosZ().y());
				if(view[i].nHits()>2)
					printf("Hit 2 position is (%lf , %lf , %lf \n) ", view[i].hitPosX().z(),view[i].hitPosY().z(),view[i].hitPosZ().z());
			}
		}
	};

	//--- Kernel for performing the seed to SC match
	//--- Still dummy and just testing out some things
	class SeedToSuperClusterMatcher {
	public:
		template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
		ALPAKA_FN_ACC void operator()(TAcc const& acc,
										reco::EleSeedDeviceCollection::View viewEleSeeds,
										int32_t sizeEleSeeds,
										reco::SuperclusterDeviceCollection::View viewSCs,
										int32_t sizeSCs,
										double vtx_x,
										double vtx_y,
										double vtx_z) const 
		{

			const int32_t thread = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0u];

			Vector3f vertex = {vtx_x,vtx_y,vtx_z}; 
			for (int32_t i : uniform_elements(acc,sizeEleSeeds)) 
			{
				Vector3f x2 = {0,0,0};
				Vector3f p2 = {0,0,0}; 

				for(int32_t j : uniform_elements(acc,sizeSCs)) 
				{
					Vector3f positionSC = {viewSCs[j].scR() * sin(viewSCs[j].scSeedTheta()) * cos(viewSCs[j].scPhi()),viewSCs[j].scR() * sin(viewSCs[j].scSeedTheta()) * cos(viewSCs[j].scPhi()),viewSCs[j].scR() * cos(viewSCs[j].scSeedTheta())};
					for (int charge : {1,-1}) 
					{
						auto newfreeTS = ftsFromVertexToPointPortable::ftsFromVertexToPoint(positionSC,vertex, viewSCs[j].scEnergy(),charge,magneticFieldParabolicPortable::magneticFieldAtPoint(positionSC));									
						Vector3f position = {newfreeTS.position(0),newfreeTS.position(1),newfreeTS.position(2)};
						Vector3f momentum = {newfreeTS.momentum(0),newfreeTS.momentum(1),newfreeTS.momentum(2)};
						
						auto transverseCurvature = [](const Vector3f& p, int charge, const float& magneticFieldZ) -> float {
   			 				return -2.99792458e-3f * (charge / sqrt(p(0) * p(0) + p(1) * p(1))) * magneticFieldZ;  
						};
						double s=0; 
						double rho = transverseCurvature(momentum,1,magneticFieldParabolicPortable::magneticFieldAtPoint(positionSC));
						bool theSolExists = false; 

						Vector3f surfPosition = {viewEleSeeds[i].surfPosX().x(),viewEleSeeds[i].surfPosY().x(),viewEleSeeds[i].surfPosZ().x()};
						Vector3f surfRotation = {viewEleSeeds[i].surfRotX().x(),viewEleSeeds[i].surfRotY().x(),viewEleSeeds[i].surfRotZ().x()};

						Propagators::helixBarrelPlaneCrossing(position,momentum,rho,Propagators::oppositeToMomentum,surfPosition,surfRotation,theSolExists,x2,p2,s);
						printf("Position : %lf and direction : %lf and path length : %lf",x2(0),p2(0),s);
					}

					printf("positionSC %lf  : %lf and direction : %lf",positionSC(0),x2(0),p2(0));
					printf(" Vertex position : (%lf, %lf, %lf)\n",vertex(0),vertex(1),vertex(2));
				}
			}
		}
	};


	//---- Kernel launch for printing the SC SoA collection
	void PixelMatchingAlgo::printEleSeeds(Queue& queue, reco::EleSeedDeviceCollection& collection) const {
		uint32_t items = 32;
		uint32_t groups = divide_up_by(collection->metadata().size(), items);
		auto workDiv = make_workdiv<Acc1D>(groups, items);
		alpaka::exec<Acc1D>(queue, workDiv, printElectronSeedSoA{}, collection.view(), collection->metadata().size());
	}

	//---- Kernel launch for printing the SC SoA collection
	void PixelMatchingAlgo::printSCs(Queue& queue, reco::SuperclusterDeviceCollection& collection) const {
		uint32_t items = 32;
		uint32_t groups = divide_up_by(collection->metadata().size(), items);
		auto workDiv = make_workdiv<Acc1D>(groups, items);
		alpaka::exec<Acc1D>(queue, workDiv, printSCSoA{}, collection.view(), collection->metadata().size());
	}

	//---- Kernel launch for SC and seed matching
	void PixelMatchingAlgo::matchSeeds(Queue& queue, reco::EleSeedDeviceCollection& collection, reco::SuperclusterDeviceCollection& collectionSCs, double vtx_X, double vtx_Y, double vtx_Z) const 
	{
		uint32_t items = 32;
		uint32_t groups = divide_up_by(collection->metadata().size(), items);
		auto workDiv = make_workdiv<Acc1D>(groups, items);
		alpaka::exec<Acc1D>(queue, workDiv, SeedToSuperClusterMatcher{},collection.view(),collection->metadata().size(),collectionSCs.view(),collectionSCs->metadata().size(),vtx_X, vtx_Y,vtx_Z);
	}

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
