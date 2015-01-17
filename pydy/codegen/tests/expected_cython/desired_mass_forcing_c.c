#include <math.h>
#include "desired_mass_forcing_c.h"

void mass_forcing(double constants[4], // constants = [k0, m0, g, c0]
                  double coordinates[1], // coordinates = [x0]
                  double speeds[1], // speeds = [v0]
                  double specified[1], // specified = [f0]
                  double mass_matrix[4], // computed
                  double forcing_vector[2]) // computed
{
    // common subexpressions
    double z_0 = speeds[0];

    // mass matrix
    mass_matrix[0] = 1;
    mass_matrix[1] = 0;
    mass_matrix[2] = 0;
    mass_matrix[3] = constants[1];

    // forcing vector
    forcing_vector[0] = z_0;
    forcing_vector[1] = -constants[3]*z_0 + constants[2]*constants[1] - constants[0]*coordinates[0] + specified[0];
}
