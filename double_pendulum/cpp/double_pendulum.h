#ifndef DOUBLE_PENDULUM_H
#define DOUBLE_PENDULUM_H
double double_pendulum_ke(const double x[4], const double params[3]);
double double_pendulum_pe(const double x[4], const double params[3]);
int double_pendulum_ode(double t, const double x[], double dxdt[], void *_params);
#endif
