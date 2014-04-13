#!/usr/bin/env python

"""This module contains templates for any files that have to be generated."""

c_template = \
"""\
#include <math.h>
#include "{header_filename}"

void mass_forcing(double constants[{constants_len}], // constants = [{constants_list}]
                  double coordinates[{coordinates_len}], // coordinates = [{coordinates_list}]
                  double speeds[{speeds_len}], // speeds = [{speeds_list}]
{specified_double}                  double mass_matrix[{mass_matrix_len}], // computed
                  double forcing_vector[{forcing_vector_len}]) // computed
{{
    // common subexpressions
    {sub_expression_block}

    // mass matrix
    {mass_matrix_block}

    // forcing vector
    {forcing_vector_block}
}}
"""

h_template = \
"""\
void mass_forcing(double constants[{constants_len}], // constants = [{constants_list}]
                  double coordinates[{coordinates_len}], // coordinates = [{coordinates_list}]
                  double speeds[{speeds_len}], // speeds = [{speeds_list}]
{specified_double}                  double mass_matrix[{mass_matrix_len}], // computed
                  double forcing_vector[{forcing_vector_len}]); // computed
"""

pyx_template = \
"""\
import numpy as np
cimport numpy as np

cdef extern from "{header_filename}":
    void mass_forcing(double* constants,
                      double* coordinates,
                      double* speeds,{cdef_specified_arg}
                      double* mass_matrix,
                      double* forcing_vector)


def mass_forcing_matrices(np.ndarray[np.double_t, ndim=1, mode='c'] constants,
                          np.ndarray[np.double_t, ndim=1, mode='c'] coordinates,
                          np.ndarray[np.double_t, ndim=1, mode='c'] speeds{def_specified_arg}):

    assert len(constants) == {constants_len}
    assert len(coordinates) == {coordinates_len}
    assert len(speeds) == {speeds_len}{specified_assert}

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] mass_matrix = np.zeros({mass_matrix_len})
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] forcing_vector = np.zeros({forcing_vector_len})

    mass_forcing(<double*> constants.data,
                 <double*> coordinates.data,
                 <double*> speeds.data,{call_specified_arg}
                 <double*> mass_matrix.data,
                 <double*> forcing_vector.data)

    return mass_matrix.reshape({forcing_vector_len}, {forcing_vector_len}), forcing_vector.reshape({forcing_vector_len}, 1)
"""
# TODO : Add a doc string to the cython function

setup_template = \
"""\
import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_module = Extension(name="{prefix}",
                       sources=["{pyx_filename}",
                                "{c_filename}"],
                       include_dirs=[numpy.get_include()])

setup(name="{prefix}",
      cmdclass = {{'build_ext': build_ext}},
      ext_modules = [ext_module])
"""
