#include "RecoEgamma/EgammaElectronAlgos/interface/TrajSeedMatcherCUDA.h"

//CUDA dependencies and libaries
#include <cuda_runtime.h>
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

//C++ libaries
#include <stdio.h>
#include <vector>
#include <cmath>

namespace TrajSeedMatcherCUDA{

//Kernel(s)
  __global__ void match(float *binArr, float &cutVal, const float etaVal, int binSize){
 

  }
  
  
  #ifdef __CUDACC__
    void wrapper(){}     
  #endif
  
}
