#include <cuda_runtime.h>
#include <stdio.h>
#include <iostream>
#include "RecoEgamma/EgammaElectronProducers/interface/TestCUDA.h"

namespace EgammaCUDA {

    __global__ void hello_world_gpu() {
        printf("Hello World from the GPU ");
        if ( blockIdx.x < 100 && threadIdx.x < 100 ) 
            printf("Hello World from the GPU at block %u, thread %u \n", blockIdx.x, threadIdx.x);
    }

    __global__ void printEcalSCkernel(unsigned int nSCs, EcalSC::EcalSCSoA* SCs) {
        size_t firstElement = threadIdx.x; 
        for (unsigned int isc = firstElement; isc < nSCs ; isc += blockDim.x){
            printf("Hello World from the GPU at block %u, thread %u \n", blockIdx.x, threadIdx.x);
            printf("Block dimention is %u \n",blockDim.x );
            printf("superClusRef->seed()->position().theta() : %lf \n",SCs->scTheta(isc) );
        }
        __syncthreads();
    }


    #ifdef __CUDACC__
        void hello_world_gpu_Wrapper() 
        {
            /* Call GPU function */
            const int n_blocks  = 1;
            const int n_threads = 32;
            dim3 grid_dim(n_blocks);
            dim3 block_dim(n_threads);
            std::cout<<" I am in ElectronNHitSeedProducerCUDA::produce "<<std::endl;
            hello_world_gpu<<<grid_dim, block_dim>>>();
            cudaDeviceSynchronize();
            cudaError_t error = cudaGetLastError();
            if (error != cudaSuccess) {
                fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
            }
        }

        void ecalScWrapper(unsigned int nSCs, EcalSC::EcalSCSoA* SCs, cudaStream_t stream)
        {
            unsigned int blockSize = 32; //Must be less that 1024 for the T4
            unsigned int gridSize  = (nSCs*blockSize-1)/blockSize; //A good practice is size*blockSize-1 / blockSize 
            printf("Number of SCs : %u \n", nSCs);
            printEcalSCkernel<<<gridSize, blockSize,0,stream>>>(nSCs, SCs);
            cudaDeviceSynchronize();
            cudaError_t error = cudaGetLastError();
            if (error != cudaSuccess) {
                fprintf(stderr, "ERROR: %s \n", cudaGetErrorString(error));
            }

            std::cout<<" Exit the device "<<std::endl;
        }

    #endif
} // End of namespace