#include "double_pendulum_ode.h"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

int double_pendulum_ode(double t, const double x[], double dxdt[], void *_params)
{
  double *params = static_cast<double *>(_params);
  double l = params[0], m = params[1], g = params[2];
  dxdt[0] = x[2];
  dxdt[1] = x[3];
  dxdt[2] = (-g*sin(x[0])*pow(sin(x[1]), 2) + 2*g*sin(x[0]) - g*sin(x[1])*cos(x[0])*cos(x[1]) + 2*l*pow(x[2], 2)*sin(x[0])*cos(x[0])*pow(cos(x[1]), 2) - l*pow(x[2], 2)*sin(x[0])*cos(x[0]) - 2*l*pow(x[2], 2)*sin(x[1])*pow(cos(x[0]), 2)*cos(x[1]) + l*pow(x[2], 2)*sin(x[1])*cos(x[1]) + l*pow(x[3], 2)*sin(x[0])*cos(x[1]) - l*pow(x[3], 2)*sin(x[1])*cos(x[0]))/(l*(pow(sin(x[0]), 2)*pow(sin(x[1]), 2) + 2*sin(x[0])*sin(x[1])*cos(x[0])*cos(x[1]) + pow(cos(x[0]), 2)*pow(cos(x[1]), 2) - 2)) ;
  dxdt[3] = (-sin(x[0])*sin(x[1])/2 - cos(x[0])*cos(x[1])/2)*(2*g*l*m*sin(x[0]) - pow(l, 2)*m*(-sin(x[0])*cos(x[1]) + sin(x[1])*cos(x[0]))*pow(x[3], 2))/(pow(l, 2)*m*(sin(x[0])*sin(x[1])/2 + cos(x[0])*cos(x[1])/2)*(sin(x[0])*sin(x[1]) + cos(x[0])*cos(x[1])) - pow(l, 2)*m) + (g*l*m*sin(x[1]) - pow(l, 2)*m*(sin(x[0])*cos(x[1]) - sin(x[1])*cos(x[0]))*pow(x[2], 2))/(pow(l, 2)*m*(sin(x[0])*sin(x[1])/2 + cos(x[0])*cos(x[1])/2)*(sin(x[0])*sin(x[1]) + cos(x[0])*cos(x[1])) - pow(l, 2)*m);
  return GSL_SUCCESS;
}
