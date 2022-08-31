#ifndef testCUDA_h
#define testCUDA_h
#include <cstddef>
#include <cstdint>

#include "DataFormats/EgammaReco/interface/EcalSCHeterogeneous.h"

namespace EgammaCUDA {
    struct egammaParameters {
        double testParam;
    };
    void hello_world_gpu_Wrapper();
    void ecalScWrapper(unsigned int nSCs, EcalSC::EcalSCSoA* SCs, cudaStream_t stream);
}
#endif