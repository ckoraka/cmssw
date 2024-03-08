#ifndef EgammaReco_EleSeedSoA_h
#define EgammaReco_EleSeedSoA_h

#include <cstdint>
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

using Vector3D = Eigen::Matrix <double,3,1>;
using Matrix = Eigen::Matrix <double,3,1>;

namespace reco {

	// SoA layout for supercluster
	GENERATE_SOA_LAYOUT(EleSeedLayout,
		SOA_SCALAR(bool, isValid),
		SOA_COLUMN(Vector3D, startingPos),
		SOA_COLUMN(Vector3D, startingMom),
		SOA_COLUMN(int32_t, id)
	)
	using EleSeedSoA = EleSeedLayout<>;
}  // namespace reco

#endif
