#ifndef REPRODUCEKANE1982_H
#define REPRODUCEKANE1982_H

struct rattleback_params {
 double a, b, c, h, M, A, B, C, D, g;
};

struct simdata { double t, x[6], delta, alpha[3], ke; };

int rattleback_ode(double t, const double x[], double dxdt[], void *params);
void rattleback_outputs(simdata *s, rattleback_params *p);

#endif
