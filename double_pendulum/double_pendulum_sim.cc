#include <fstream>
#include <gsl/gsl_odeiv2.h>
#include "double_pendulum_ode.h"

// A Data structure to store time, state, potential and kinetic energy at each
// time step
struct simdata {
  double t;
  double x[4];
  double pe;
  double ke;
};

int main(int argc, char *argv[])
{
  using namespace std;
  // Simulation conditions
  double params[] = {1.0, 9.81, 1.0};       // m, g, l
  simdata s = {0.0, {0.1, 0.2, 0.0, 0.0}};  // initial time and state
  s.ke = double_pendulum_ke(s.x, params);
  s.pe = double_pendulum_pe(s.x, params);
  const double tf = 5.0;                    // final time
  const int N = 501;                        // number of points
  
  // GSL setup code
  gsl_odeiv2_system sys = {double_pendulum_ode, NULL, 4, params};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, 
                           gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

  // Open a file for writing
  ofstream f("datafile.dat", ios::binary | ios::out);

  // Simulation driver loop
  f.write((char *) &s, sizeof(simdata));
  for (int i = 1; i <= N; ++i) {
    double ti = i * tf / N;
    gsl_odeiv2_driver_apply(d, &(s.t), ti, s.x);
    s.ke = double_pendulum_ke(s.x, params);
    s.pe = double_pendulum_pe(s.x, params);
    f.write((char *) &s, sizeof(simdata));
  }

  gsl_odeiv2_driver_free(d);
  f.close();
  return 0;
}
