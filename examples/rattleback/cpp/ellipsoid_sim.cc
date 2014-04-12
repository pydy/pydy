#include <fstream>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include "ellipsoid_no_slip.h"

int main() {
  using namespace Mechanics::Rattleback::Ellipsoid::NoSlip;
  parameters p;
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
  p.s = 0.0001;             // sigma, viscous air damping

  // Initial time and state
  simdata s = {0.0,              // time
              {0.0,              // Yaw (ignorable)
               0.5*M_PI/180.0,   // Roll
               0.5*M_PI/180.0,   // Pitch
               0.0, 0.0,         // x, y of contact (ignorable)
               0.0,              // u0
               0.0,              // u1
               -5.0} };          // u2  (spin)

  // GSL setup code
  gsl_odeiv2_system sys = {ode, jacobian, 8, &p};
  gsl_odeiv2_driver * d =  gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0); 

  // Open a file for writing
  std::ofstream f("simulation.dat", std::ios::binary | std::ios::out);

  // Simulation
  const double tf = 60.0;                // final time
  const int N = 20001;                   // number of points
  outputs(&s, &p);           // compute initial output quantities
  f.write((char *) &s, sizeof(simdata));// Write initial time data
  for (int i = 1; i <= N; ++i) {
    double ti = i * tf / N;
    gsl_odeiv2_driver_apply(d, &(s.t), ti, s.x);  // integrate the ODE's
    outputs(&s, &p);                   // compute output quantities
    f.write((char *) &s, sizeof(simdata));        // write to file
  }

  gsl_odeiv2_driver_free(d);                      // free resources
  f.close();                                      // close file
  return 0;
}
