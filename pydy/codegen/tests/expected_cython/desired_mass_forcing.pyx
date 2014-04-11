import numpy as np
cimport numpy as np

cdef extern from "desired_mass_forcing_c.h":
    void mass_forcing(double* constants,
                      double* coordinates,
                      double* speeds,
                      double* specified,
                      double* mass_matrix,
                      double* forcing_vector)


def mass_forcing_matrices(np.ndarray[np.double_t, ndim=1, mode='c'] constants,
                          np.ndarray[np.double_t, ndim=1, mode='c'] coordinates,
                          np.ndarray[np.double_t, ndim=1, mode='c'] speeds,
                          np.ndarray[np.double_t, ndim=1, mode='c'] specified):

    assert len(constants) == 4
    assert len(coordinates) == 1
    assert len(speeds) == 1
    assert len(specified) == 1


    cdef np.ndarray[np.double_t, ndim=1, mode='c'] mass_matrix = np.zeros(4)
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] forcing_vector = np.zeros(2)

    mass_forcing(<double*> constants.data,
                 <double*> coordinates.data,
                 <double*> speeds.data,
                 <double*> specified.data,
                 <double*> mass_matrix.data,
                 <double*> forcing_vector.data)

    return mass_matrix.reshape(2, 2), forcing_vector.reshape(2, 1)
