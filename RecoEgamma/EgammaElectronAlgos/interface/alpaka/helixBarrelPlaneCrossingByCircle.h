/**
 Description: Function to propagate from a point to a plane on the GPU
*/

#ifndef RecoEgamma_EgammaElectronAlgos_interface_alpaka_helixBarrelPlaneCrossing_h
#define RecoEgamma_EgammaElectronAlgos_interface_alpaka_helixBarrelPlaneCrossing_h

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

#include "DataFormats/EgammaReco/interface/Plane.h"

using Vector3D = Eigen::Matrix<float, 3, 1>;

//enum PropagationDirection { alongMomentum, oppositeMomentum, anyDirection };

namespace ALPAKA_ACCELERATOR_NAMESPACE {

	namespace Propagators {

		ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE void helixBarrelPlaneCrossing(
		const Vector3D& startingPos,
		const Vector3D& startingDir,
		double rho,
		PropagationDirection propDir,
		Vector3D& surfPosition,
		Vector3D& surfRotation,            
		bool& theSolExists,
		Vector3D& position,
		Vector3D& direction,
		double& s) 
		{

		const double sraightLineCutoff = 1.e-7;

		double pt = startingDir.head(2).norm();
		double ipabs = 1. / startingDir.norm();
		double sinTheta = pt * ipabs;
		double cosTheta = startingDir(2) * ipabs;

		// === TO FIX SHOULD ADD STRAIGHT LINE CASE === //
		//bool useStraightLine;
		if (fabs(rho) < sraightLineCutoff && fabs(rho) * startingPos.head(2).norm() < sraightLineCutoff) {
			//useStraightLine = true;
			// Add implementation
			//StraightLinePlaneCrossing slc(startingPos, startingDir, propDir);
			//auto pl = slc.pathLength(plane);
			//if (pl.first) {
			//	theSolExists = true;
			//	s = pl.second;
			//	position = slc.position(s);
			//	direction = startingDir;
			//} else {
				theSolExists = false;
			//}
			//return; // all needed data members have been set
		}

		double o = 1. / (pt * rho);
		double theXCenter = startingPos(0) - startingDir(1) * o;
		double theYCenter = startingPos(1) + startingDir(0) * o;

		// This is default when there curvature is non zero

		PlanePortable::Plane<Vector3D> plane{surfPosition,surfRotation};
		Vector3D n = plane.normalVector();
		double distToPlane = -plane.localZ(startingPos);
		double nx = n(0);
		double ny = n(1);
		double distCx = startingPos(0) - theXCenter;
		double distCy = startingPos(1) - theYCenter;

		double nfac, dfac;
		double A, B, C;
		bool solveForX;

		if (fabs(nx) > fabs(ny)) {
			solveForX = false;
			nfac = ny / nx;
			dfac = distToPlane / nx;
			B = distCy - nfac * distCx; // only part of B
			C = (2. * distCx + dfac) * dfac;
		} else {
			solveForX = true;
			nfac = nx / ny;
			dfac = distToPlane / ny;
			B = distCx - nfac * distCy; // only part of B
			C = (2. * distCy + dfac) * dfac;
		}

		B -= nfac * dfac;
		B *= 2; // the rest of B
		A = 1. + nfac * nfac;

		QuadEquationSolver eq(A, B, C);
		if (!eq.hasSolution) {
			theSolExists = false;
			return;
		}

		Vector3D d1, d2;
		if (solveForX) {
			d1 = Vector3D(eq.first, dfac - nfac * eq.first, 0.0);
			d2 = Vector3D(eq.second, dfac - nfac * eq.second, 0.0);
		} else {
			d1 = Vector3D(dfac - nfac * eq.first, eq.first, 0.0);
			d2 = Vector3D(dfac - nfac * eq.second, eq.second, 0.0);
		}

		Vector3D theD;
		int theActualDir;
		theD = chooseSolution(d1, d2, startingPos, startingDir, propDir, theActualDir, theSolExists);
		if (!theSolExists)
			return;

		auto dMag = theD.norm();
		double sinAlpha = 0.5 * dMag * rho;
		if (std::abs(sinAlpha) > 1.0)
			 sinAlpha = std::copysign(1.0, sinAlpha);
		// Path length
		s = theActualDir * 2.0 * std::asin(sinAlpha) / (rho * sinTheta);

		// Position
		position = Vector3D(startingPos(0) + theD(0), startingPos(1) + theD(1), startingPos(2) + s * cosTheta);

		// Direction
		double sinPhi, cosPhi;
    	double tmp = 0.5 * dMag * rho;
		if (s < 0) 
			tmp = -tmp;
 	    sinPhi = 1. - (tmp * tmp);
		if (sinPhi < 0)
			sinPhi = 0.;
		sinPhi = 2.0 * tmp * sqrt(sinPhi);
		cosPhi = 1.0 - 2.0 * (tmp * tmp);

		direction = Vector3D(startingDir(0) * cosPhi - startingDir(1) * sinPhi,
							  startingDir(0) * sinPhi + startingDir(1) * cosPhi,
							  startingDir(2));
		}

		ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE Vector3D chooseSolution(
			const Vector3D& d1,
			const Vector3D& d2,
			const Vector3D& startingPos,
			const Vector3D& startingDir,
			PropagationDirection propDir,
			int& theActualDir,
			bool& theSolExists) {

			Vector3D theD;

			auto momProj1 = startingDir(0) * d1(0) + startingDir(1) * d1(1);
			auto momProj2 = startingDir(0) * d2(0) + startingDir(1) * d2(1);

			if (propDir == anyDirection) {
				theSolExists = true;
				if (d1.squaredNorm() < d2.squaredNorm()) {
					theD = d1;
					theActualDir = (momProj1 > 0) ? 1 : -1;
				} else {
					theD = d2;
					theActualDir = (momProj2 > 0) ? 1 : -1;
				}
			} else {
				double propSign = (propDir == alongMomentum) ? 1 : -1;
				if (momProj1 * momProj2 < 0) {
					theSolExists = true;
					theD = (momProj1 * propSign > 0) ? d1 : d2;
					theActualDir = propSign;
				} else if (momProj1 * propSign > 0) {
					theSolExists = true;
					theD = (d1.squaredNorm() < d2.squaredNorm()) ? d1 : d2;
					theActualDir = propSign;
				} else {
					theSolExists = false;
				}
			}

			return theD;
		}

	} // namespace Propagators

} // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
