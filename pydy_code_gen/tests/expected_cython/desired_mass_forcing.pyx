import numpy as np
cimport numpy as np

cdef extern from "desired_mass_forcing.h":
    void mass_forcing(double constants[4],
                      double coordinates[1],
                      double speeds[1],
                      double specified[1],
                      double mass_matrix[4],
                      double forcing_vector[2])


def mass_forcing_matrices(np.ndarray[np.double_t, ndim=1] constants,
                          np.ndarray[np.double_t, ndim=1] coordinates,
                          np.ndarray[np.double_t, ndim=1] speeds,
                          np.ndarray[np.double_t, ndim=1] specified):

    cdef np.ndarray[np.double_t, ndim=1] mass_matrix = np.zeros(4)
    cdef np.ndarray[np.double_t, ndim=1] forcing_vector = np.zeros(2)

    mass_forcing(<double*> constants.data,
                 <double*> coordinates.data,
                 <double*> speeds.data,
                 <double*> specified.data,
                 <double*> mass_matrix.data,
                 <double*> forcing_vector.data)

    return mass_matrix.reshape(4, 1), forcing_vector.reshape(2, 1)
