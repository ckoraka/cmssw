#include <iostream>
#include <cmath>
#include <cfloat>
#include <vdt/vdtMath.h>

#include "DataFormats/GeometryVector/interface/Basic3DVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "FWCore/Utilities/interface/isFinite.h"

typedef GlobalPoint PositionType;
typedef GlobalVector DirectionType;

//Double precision vectors
typedef Basic3DVector<double> PositionTypeDouble;
typedef Basic3DVector<double> DirectionTypeDouble;

//Functions and structures
PositionTypeDouble updatePos(double s, double theCachedS, double theCachedDPhi, double theCachedCDPhi, double theCachedSDPhi, double sinPhi0, double cosPhi0, double rho, double cosTheta, double sinTheta, double sinThetaI, double x0, double y0, double z0){
  double xNew, yNew, zNew;
  if (s != theCachedS){
    theCachedS = s;
    theCachedDPhi = theCachedS * rho * sinTheta;
    vdt::fast_sincos(theCachedDPhi, theCachedSDPhi, theCachedCDPhi);
  }
  
  if (std::abs(theCachedDPhi) > 1.e-4f){
    double o = 1./rho;
    xNew = x0 + (-sinPhi0 * (1. - theCachedCDPhi) + cosPhi0 * theCachedSDPhi) * o;
    yNew = y0 + (cosPhi0 * (1. - theCachedCDPhi) + sinPhi0 * theCachedSDPhi) * o;
    zNew = z0 + theCachedS * cosTheta; 
  }else{
    //Use second order
    double st = s/sinThetaI;
    xNew = x0 + (cosPhi0 - (st * 0.5 * rho) * sinPhi0) * st;
    yNew = y0 + (sinPhi0 + (st * 0.5 * rho) * cosPhi0) * st;
    zNew = z0 + st * cosTheta * sinThetaI;
  }

  return PositionTypeDouble(xNew, yNew, zNew);
}

std::pair<bool, double> SolutionByDirection(const double dS1, const double dS2, int thePropDir){
  bool valid = false;
  double path = 0;
  if (thePropDir == 0) { //How do I deal with thePropDir??
    valid = true;
    path = fabs(dS1) < fabs(dS2) ? dS1 : dS2;
  } else {
    // use same logic for both directions (invert if necessary)
    double propSign = thePropDir == 1 ? 1 : -1;
    double s1(propSign * dS1);
    double s2(propSign * dS2);
    // sort
    if (s1 > s2)
      std::swap(s1, s2);
    // choose solution (if any with positive sign)
    if ((s1 < 0) & (s2 >= 0)) {
      // First solution in backward direction: choose second one.
      valid = true;
      path = propSign * s2;
    } else if (s1 >= 0) {
      // First solution in forward direction: choose it (s2 is further away!).
      valid = true;
      path = propSign * s1;
    }
  }
  if (edm::isNotFinite(path))
    valid = false;
  return std::pair<bool, double>(valid, path);
}

std::pair<bool, double> PathLenghtFirstPass(const PositionType& point, const DirectionType& direction, const float curvature, const PositionType& planePos, const GlobalVector& normalVec, int const propDir){
  double theRho = curvature;
  
  //Break up the normal vector
  double nPx = normalVec.x();
  double nPy = normalVec.y();
  double nPz = normalVec.z();
  double cP = normalVec.dot(point - planePos);
   
  // Components of direction vector (with correct normalisation)
  double px = direction.x();
  double py = direction.y();
  double pz = direction.z();
  double pt2 = px * px + py * py;
  double p2 = pt2 + pz * pz;
  double pI = 1. / sqrt(p2);
  double ptI = 1. / sqrt(pt2);
  double theCosPhi0 = px * ptI;
  double theSinPhi0 = py * ptI;
  double theCosTheta = pz * pI;
  double theSinThetaI = p2 * pI * ptI;  //  (1/(pt/p)) = p/pt = p*ptI and p = p2/p = p2*pI

  //Coefficients of 2nd order equation to obtain interseciton point (no curvature-related factors)
  double ceq1 = theRho * (nPx * theSinPhi0 - nPy * theCosPhi0);
  double ceq2 = nPx * theCosPhi0 + nPy * theSinPhi0 + nPz * theCosTheta * theSinThetaI;
  double ceq3 = cP;
  
  // Check for degeneration to linear equation (zero
  //   curvature, forward plane or direction perp. to plane) 
  double dS1, dS2;
  if (std::abs(ceq1) > FLT_MIN) {
    double deq1 = ceq2 * ceq2;
    double deq2 = ceq1 * ceq3;
    if (std::abs(deq1) < FLT_MIN || std::abs(deq2 / deq1) > 1.e-6) {
      
      // Standard solution for quadratic equations
      double deq = deq1 + 2 * deq2;
      if (deq < 0.)
        return std::pair<bool, double>(false, 0);
      double ceq = ceq2 + std::copysign(std::sqrt(deq), ceq2);
      dS1 = (ceq / ceq1) * theSinThetaI;
      dS2 = -2. * (ceq3 / ceq) * theSinThetaI;
    } else {
      //
      // Solution by expansion of sqrt(1+deq)
      //
      double ceq = (ceq2 / ceq1) * theSinThetaI;
      double deq = deq2 / deq1;
      deq *= (1 - 0.5 * deq);
      dS1 = -ceq * deq;
      dS2 = ceq * (2 + deq);
    }
  } else {
    
    // Special case: linear equation
    dS1 = dS2 = -(ceq3 / ceq2) * theSinThetaI;
  }

  // Choose and return solution
  return SolutionByDirection(dS1, dS2, propDir);
}


//Propagator
void helixToArbitraryPlanePropagator(const PositionType& point,
                                     const DirectionType& direction,
                                     const double curvature,
                                     const PositionType& planePos,
                                     const GlobalVector& normalVec,
                                     double& theS,
                                     bool& solution,
                                     GlobalPoint& x, GlobalVector& p){
  double theX0 = point.x();
  double theY0 = point.y();
  double theZ0 = point.z();
  double theRho = curvature;
  double theCachedS = 0.;
  double theCachedDPhi = 0.;
  double theCachedSDPhi = 0.;
  double theCachedCDPhi = 1.;

  //Normalized direction vectors
  double px = direction.x();
  double py = direction.y();
  double pz = direction.z();
  double pt2 = px * px + py * py;
  double p2 = pt2 + pz * pz;
  double pI = 1. / sqrt(p2);
  double ptI = 1. / sqrt(pt2);
  double theCosPhi0 = px * ptI;
  double theSinPhi0 = py * ptI;
  double theCosTheta = pz * pI;
  double theSinTheta = pt2 * ptI * pI;
  double theSinThetaI = p2 * pI * ptI;  //  (1/(pt/p)) = p/pt = p*ptI and p = p2/p = p2*pI

  //Path Length and solution/notFailure information
  bool notFail = solution;

  //Convergance control and constraints
  const float numPrecision = 5.e-7f;
  const float theMaxDistToPlane = 1.e-4f;

  //Maximum interation limit
  //constexpr int maxIterations = 20;

  //Distance to plane
  float maxNumDZ = numPrecision*planePos.mag();
  float safeMaxDist = (theMaxDistToPlane > maxNumDZ ? theMaxDistToPlane : maxNumDZ);
  float dz = normalVec.dot(point - planePos); //plane.localZ equivelent
  int propDir = 1; //set prop dir along direction of momenta
  
  if (std::abs(dz) < safeMaxDist){
    std::cout << "0 propagation" << std::endl;
    solution = true;
    theS = 0;
    return;
  }
 
  double dSTotal;
  //Obtain inital pass for path lenght using second order object
  std::tie(notFail, dSTotal) = PathLenghtFirstPass(point, direction, curvature, planePos, normalVec, propDir);
  if (!notFail){
    std::cout << "No Solution" << std::endl;
    solution = notFail;
    theS = 0;
    return;
  }
  
  //Update direction
  PositionTypeDouble newpos = updatePos(dSTotal, theCachedS, theCachedDPhi, theCachedCDPhi, theCachedSDPhi, theSinPhi0, theCosPhi0, theRho, theCosTheta, theSinTheta, theSinThetaI, theX0, theY0, theZ0); 
  
  auto newDir = dSTotal >=0 ? 1 : -1; //Forward propogaton if dSTotal >= 0
  if (newDir != propDir){
     std::cout << "No Solution" << std::endl;
     solution = false;
     theS = 0;
     return; 
  }

  //Iterations count and total path length
  //auto iteration = maxIterations;
  DirectionType pnew;
  while (theMaxDistToPlane < std::abs(normalVec.dot(PositionType(newpos.x(), newpos.y(), newpos.z()) - planePos))){
    double dph = dSTotal*theRho/theSinTheta;
    /*if (iteration--){
      std::cout << "Failed to reach solution within " << maxIterations << " interations" << std::endl;
      return;
    }
    iteration--;*/
    
    //Obtain new direction
    double pxNew, pyNew, pzNew;
    if (dSTotal != theCachedS){
      theCachedS = dSTotal;
      theCachedDPhi = theCachedS * theRho * theSinTheta;
      vdt::fast_sincos(theCachedDPhi, theCachedSDPhi, theCachedCDPhi);
    }
  
    if (std::abs(theCachedDPhi) > 1.e-4f){
      double o = 1./theRho;
      pxNew = theCosPhi0 * theCachedCDPhi - theSinPhi0 * theCachedSDPhi;
      pyNew = theSinPhi0 * theCachedCDPhi + theCosPhi0 * theCachedSDPhi;
      pzNew = theCosTheta / theSinTheta; 
    }else{
      //Use second order
      double dph = dSTotal* theRho/ theSinThetaI;
      pxNew = theCosPhi0 - (theSinPhi0 + 0.5 * dph * theCosPhi0) * dph; 
      pyNew = theSinPhi0 + (theCosPhi0 - 0.5 * dph * theSinPhi0) * dph;
      pzNew = theCosTheta * theSinThetaI;
    }
    
    pnew = DirectionType(pxNew, pyNew, pzNew); 
    auto deltaS2 = PathLenghtFirstPass(PositionType(newpos.x(), newpos.y(), newpos.z()), DirectionType(pxNew, pyNew, pzNew), theRho, planePos, normalVec, propDir);
    if (!deltaS2.first){
      std::cout << "No Solution" << std::endl;
      solution = deltaS2.first;
      theS = 0;
      return;
    }
    
    //Obtain total path length
    dSTotal += deltaS2.second;
    newDir = dSTotal >= 0 ? 1 : -1; 
    if (newDir != propDir){
      std::cout << "No Solution" << std::endl;
      solution = false;
      theS = 0;
      return;
    }
   
   
   newpos = updatePos(dSTotal, theCachedS, theCachedDPhi, theCachedCDPhi, theCachedSDPhi, theSinPhi0, theCosPhi0, theRho, theCosTheta, theSinTheta, theSinThetaI, theX0, theY0, theZ0);
  
  } 

  //Update path length momenta and position
  theS = dSTotal;
  x = PositionType(newpos.x(), newpos.y(), newpos.z());
  p = pnew;
  
}
