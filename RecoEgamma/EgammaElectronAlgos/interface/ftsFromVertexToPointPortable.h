#ifndef RecoEgamma_EgammaElectronAlgos_ftsFromVertexToPointPortable_h
#define RecoEgamma_EgammaElectronAlgos_ftsFromVertexToPointPortable_h

#include <cmath>
#include <Eigen/Core>

namespace ftsFromVertexToPointPortable {

    // FreeTrajectoryState template structure
    template <typename Vec3>
    struct FreeTrajectoryState {
        Vec3 momentum;  // 3D momentum vector
        Vec3 position;  // 3D position vector
        int charge;     // Particle charge

        // Constructor
        FreeTrajectoryState(Vec3 p, Vec3 pos, int q) 
            : momentum(p), position(pos), charge(q) {}
    };

    // Function to calculate the FreeTrajectoryState from vertex to point
    template <typename Vec3>
    FreeTrajectoryState<Vec3> ftsFromVertexToPoint(
        const Vec3& xmeas,    // Measured point
        const Vec3& xvert,    // Vertex point
        float momentum,       // Magnitude of momentum
        int charge,           // Charge of the particle
        float BInTesla        // Magnetic field (in Tesla)
    ) {
        // Calculate the difference between measurement and vertex positions
        Vec3 xdiff = xmeas - xvert;

        // Normalize xdiff and scale by momentum to get the momentum vector
        Vec3 mom = momentum * (xdiff / sqrt(xdiff(0) * xdiff(0) + xdiff(1) * xdiff(1) + xdiff(2) * xdiff(2)));

        // Transverse momentum (perpendicular to the z-axis)
        float pt = sqrt(mom(0) * mom(0) + mom(1) * mom(1));
        float pz = mom(2);

        float pxOld = mom(0);
        float pyOld = mom(1);

        // Calculate the curvature (assuming charge is either +1 or -1)
        float curv = (BInTesla * 0.29979f * 0.01f) / pt;

        // Calculate the sine and cosine of the rotation angle
        float sa = 0.5f * sqrt(xdiff(0) * xdiff(0) + xdiff(1) * xdiff(1)) * curv * float(charge);
        float ca = sqrt(1.f - sa * sa);

        // Rotate momentum vector in the xy-plane
        float pxNew = ca * pxOld + sa * pyOld;
        float pyNew = -sa * pxOld + ca * pyOld;
        Vec3 pNew(pxNew, pyNew, pz);

        return FreeTrajectoryState<Vec3>(pNew, xmeas, charge);
    }

}  // namespace ftsFromVertexToPointPortable

#endif  // RecoEgamma_EgammaElectronAlgos_ftsFromVertexToPointPortable_h
