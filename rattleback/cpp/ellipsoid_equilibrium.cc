#include <cmath>
#include "ellipsoid_equilibrium.h"
#include "vtkObjectFactory.h"
namespace Mechanics { namespace Rattleback { namespace Ellipsoid {
namespace NoSlip {

vtkStandardNewMacro(EllipsoidEquilibrium);

// Construct rattleback with parameters from Kane
EllipsoidEquilibrium::EllipsoidEquilibrium()
{
  p.a = 0.2;                // x semi-diameter
  p.b = 0.03;               // y semi-diameter
  p.c = 0.02;               // z semi-diameter 
  p.d = p.e = 0.0;          // COM offset in x, y directions
  p.f = 0.01;               // COM offset in z direction
  p.m = 1.0;                // mass
  p.g = 9.81;               // graviational constant
  p.Ixx =  0.0002;          // xx moment of inertia
  p.Iyy =  0.0016;          // yy moment of inertia
  p.Izz =  0.0017;          // zz moment of inertia
  p.Ixy = -0.00002;         // xy product of inertia
  p.Iyz = p.Ixz = 0.0;      // yz, xz products of inertia
  p.s = 0.0;             // sigma, viscous air damping
}

double EllipsoidEquilibrium::EvaluateFunction(double x[3])
{
  return 0.0;
}

void EllipsoidEquilibrium::EvaluateGradient(double x[3], double g[3])
{

}

void EllipsoidEquilibrium::SetParameters(parameters &params)
{
  p = params;
}
}}}}
