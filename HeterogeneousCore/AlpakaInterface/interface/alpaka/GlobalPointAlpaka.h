#ifndef HeterogeneousCore_AlpakaInterface_interface_alpaka_GlobalPointAlpaka_h
#define HeterogeneousCore_AlpakaInterface_interface_alpaka_GlobalPointAlpaka_h

#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include <type_traits>

using namespace alpaka_common;

class GlobalPointAlpaka {
   public:
       GlobalPointAlpaka(Vec1D point);
};

#endif  // HeterogeneousCore_AlpakaInterface_interface_alpaka_GlobalPointAlpaka_h