#ifndef RATTLEBACK_ELLIPSOID_NO_SLIP_H
#define RATTLEBACK_ELLIPSOID_NO_SLIP_H

namespace Mechanics { namespace Rattleback { namespace Ellipsoid {
namespace NoSlip {
  /// Constant parameters that define the Rattleback geometry, mass,
  /// mass distribution, damping coefficient and gravitational constant
  struct parameters {
   double a, b, c, d, e, f, m, Ixx, Iyy, Izz, Ixy, Iyz, Ixz, g, s;
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

  /// Dynamic equilibrium equations, only depend on roll, pitch, and yaw rate
  void ode_eq(const double roll_pitch_yawrate[], double f_dyn_eq[], void *params);
}
}}}
#endif
