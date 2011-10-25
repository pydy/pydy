#include <fstream>
#include <gsl/gsl_odeiv2.h>
#include "double_pendulum.h"

struct simdata { double t, x[4], pe, ke; };       // t, q1, q2, u1, u2, pe, ke

int main(int argc, char *argv[]) {
  double params[] = {1.0, 9.81, 1.0};             // m, g, l
  simdata s = {0.0, {0.1, 0.2, 0.0, 0.0}};        // initial time and state
  s.ke = double_pendulum_ke(s.x, params);         // initial kinetic energy
  s.pe = double_pendulum_pe(s.x, params);         // intial potential energy
  const double tf = 5.0;                          // final time
  const int N = 501;                              // number of points
  
  // GSL setup code
  gsl_odeiv2_system sys = {double_pendulum_ode, NULL, 4, params};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

  // Open a file for writing
  std::ofstream f("datafile.dat", std::ios::binary | std::ios::out);

  // Simulation loop
  f.write((char *) &s, sizeof(simdata));          // Write initial time data
  for (int i = 1; i <= N; ++i) {
    double ti = i * tf / N;
    gsl_odeiv2_driver_apply(d, &(s.t), ti, s.x);  // integrate the ODE's
    s.ke = double_pendulum_ke(s.x, params);       // compute kinetic energy
    s.pe = double_pendulum_pe(s.x, params);       // compute potential energy
    f.write((char *) &s, sizeof(simdata));        // write to file
  }

  gsl_odeiv2_driver_free(d);                      // free resources
  f.close();                                      // close file
  return 0;
}
