#ifndef DataFormats_EgammaReco_interface_Plane_h
#define DataFormats_EgammaReco_interface_Plane_h

#include <Eigen/Dense>
#include <cmath>

using Vector3D = Eigen::Matrix<float, 3, 1>;

namespace PlanePortable {

    template <typename Vector3D>
    struct Plane {
        Vector3D position;
        Vector3D rotation; // This is rotation().z() 

        // Returns the normal vector of the plane
        constexpr Vector3D normalVector() const {
            return rotation;
        }

        // Fast access to distance from plane for a point
        constexpr float localZ(const Vector3D& vp) const {
            return normalVector().dot(vp - position);
        }

        // Clamped distance from plane for a point
        constexpr float localZclamped(const Vector3D& vp) const {
            auto d = localZ(vp);
            return std::abs(d) > 1e-7 ? d : 0; 
        }

        // Fast access to distance from plane for a vector
        constexpr float distanceFromPlaneVector(const Vector3D& gv) const {
            return normalVector().dot(gv);
        }

    };

}  // namespace PlanePortable

#endif // DataFormats_EgammaReco_interface_Plane_h

