void mass_forcing(double constants[4], // constants = [c0, g, k0, m0]
                  double coordinates[1], // coordinates = [x0]
                  double speeds[1], // speeds = [v0]
                  double specified[1], // specified = [f0]
                  double mass_matrix[4], // computed
                  double forcing_vector[2]); // computed
