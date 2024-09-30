#ifndef DataFormats_EgammaReco_interface_Plane_h
#define DataFormats_EgammaReco_interface_Plane_h

#include <Eigen/Dense>
#include <cmath>

using Vector3f = Eigen::Matrix<double, 3, 1>;

namespace PlanePortable {

    template <typename Vector3D>
    struct Plane {
        Vector3f position;
        Vector3f rotation; // This is rotation().z() 

        // Returns the normal vector of the plane
        Vector3f normalVector() const {
            return rotation;
        }

        // Fast access to distance from plane for a point
        float localZ(const Vector3f& vp) const {
            return normalVector().dot(vp - position);
        }

        // Clamped distance from plane for a point
        float localZclamped(const Vector3f& vp) const {
            auto d = localZ(vp);
            return std::abs(d) > 1e-7 ? d : 0; 
        }

        // Fast access to distance from plane for a vector
        float distanceFromPlaneVector(const Vector3f& gv) const {
            return normalVector().dot(gv);
        }

    };

}  // namespace PlanePortable

#endif // DataFormats_EgammaReco_interface_Plane_h

