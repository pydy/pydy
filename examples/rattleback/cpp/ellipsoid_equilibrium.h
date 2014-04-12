#ifndef ELLIPSOID_EQUILIBRIUM_H
#define ELLIPSOID_EQUILIBRIUM_H

#include "vtkImplicitFunction.h"
#include "ellipsoid_no_slip.h"
namespace Mechanics { namespace Rattleback { namespace Ellipsoid {
namespace NoSlip {
class EllipsoidEquilibrium : public vtkImplicitFunction
{
public:
  vtkTypeMacro(EllipsoidEquilibrium, vtkImplicitFunction);
  // Construct quadric with all coefficients = 1.
  static EllipsoidEquilibrium *New();

  // Evaluate quadric equation.
  double EvaluateFunction(double x[3]);

  // Evaluate the gradient to the quadric equation.
  void EvaluateGradient(double x[3], double g[3]);

  // Set the rattleback parameters
  void SetParameters(parameters &params);

protected:
  EllipsoidEquilibrium();
  ~EllipsoidEquilibrium() {};
  parameters p;

private:
  EllipsoidEquilibrium(const EllipsoidEquilibrium&);  // Not implemented.
  void operator=(const EllipsoidEquilibrium&);  // Not implemented.
};
}}}}
#endif
