#ifndef EgammaReco_EleSeedSoA_h
#define EgammaReco_EleSeedSoA_h

#include <cstdint>
#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

using Vector3D = Eigen::Matrix <float,3,1>;
using Matrix3D = Eigen::Matrix <float,3,3>;

namespace reco {
	GENERATE_SOA_LAYOUT(EleSeedLayout,
		SOA_SCALAR(int32_t, nHits),
		SOA_COLUMN(Matrix3D, hitPos),
		SOA_COLUMN(Matrix3D, surfPos),
		SOA_COLUMN(Matrix3D, surfRot),
		SOA_COLUMN(int32_t, detectorID),
		SOA_COLUMN(bool, isValid)
	)
	using EleSeedSoA = EleSeedLayout<>;
}  // namespace reco

#endif
