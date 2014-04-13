void mass_forcing(double constants[4], // constants = [m, k, c, g]
                  double coordinates[1], // coordinates = [x]
                  double speeds[1], // speeds = [v]
                  double specified[1], // specified = [F]
                  double mass_matrix[4], // computed
                  double forcing_vector[2]); // computed
