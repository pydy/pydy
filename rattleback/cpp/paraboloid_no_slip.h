#ifndef PARABOLOID_NO_SLIP_H
#define PARABOLOID_NO_SLIP_H

namespace Mechanics { namespace Rattleback { namespace Paraboloid {
namespace NoSlip {
  /// Constant parameters that define the Rattleback geometry, mass,
  /// mass distribution, damping coefficient and gravitational constant
  struct parameters {
   double a, b, c, d, e, m, Ixx, Iyy, Izz, Ixy, Iyz, Ixz, g, s;
  };

  /// Quantities to be calculated and saved every time step
  /// time, state, state derivatives, contact forces, kinetic and potential
  /// energy, "wobble angle"
  struct simdata { double t, x[8], dxdt[8], CF[3], ke, pe, te, delta; };

  /// Kinematic and dynamic equations of motion in first order form; given t,
  /// x, and parameters, compute dxdt
  int ode(double t, const double x[], double dxdt[], void *params);

  /// Jacobian matrix; given t, x, and parameters, compute dfdx and dfdt
  int jacobian(double t, const double x[],
               double dfdx[], double dfdt[], void *params);

  /// Output quantities; given a simdata struct with the t and x fields
  /// populated, and a set rattleback parameters, compute all output quantities
  void outputs(simdata *s, parameters *p);
}
}}}
#endif
