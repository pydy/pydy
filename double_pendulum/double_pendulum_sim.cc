#include <iostream>
#include <gsl/gsl_odeiv2.h>
#include "double_pendulum_ode.h"

int main(int argc, char *argv[])
{
  // GSL setup code
  gsl_odeiv2_system sys = {double_pendulum_ode, NULL, 4, params};
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

  // Simulation conditions
  double params[] = {1.0, 1.0, 1.0};  // l, m, g
  double x[] = {0.1, 0.2, 0.0, 0.0};  // initial conditions
  double ti = 0.0;                    // initial time
  const double tf = 10.0;             // final time
  const int N = 100;                  // number of points

  // Simulation driver loop
  for (int i = 1; i <= 100; ++i) {
    double ti = i * tf / N;
    gsl_odeiv2_driver_apply(d, &t, ti, x);
    printf("%.5e %.5e %.5e %.5e %.5e\n", t, x[0], x[1], x[2], x[3]);
  }

  gsl_odeiv2_driver_free(d);
  return 0;
}
