#!/usr/bin/env python

# standard library
import re

# external libraries
import numpy as np
from sympy import Dummy, lambdify, numbered_symbols, cse
from sympy.printing import ccode


def generate_mass_forcing_cython_code(filename_prefix, mass_matrix,
                                      forcing_vector, coordinates, speeds,
                                      specified, constants,
                                      time_variable='t'):

    # TODO : Maybe, allow expressions to be passed in for the specified
    # quantities for computation inside the c file.

    c_filename = filename_prefix + '_c.c'
    header_filename = filename_prefix + '_c.h'
    pyx_filename = filename_prefix + '.pyx'
    setup_py_filename = filename_prefix + '_setup.py'

    c_template = \
"""\
#include <math.h>
#include "{}"

void mass_forcing(double constants[{}],
                  double coordinates[{}],
                  double speeds[{}],
                  double specified[{}],
                  double mass_matrix[{}], // computed
                  double forcing_vector[{}]) // computed
{{
    // constants = [{}]
    // coordinates = [{}]
    // speeds = [{}]
    // specified = [{}]

    // common subexpressions
    {}

    // mass matrix
    {}

    // forcing vector
    {}
}}
"""

    h_template = \
"""\
void mass_forcing(double constants[{}], // constants = [{}]
                  double coordinates[{}], // coordinates = [{}]
                  double speeds[{}], // speeds = [{}]
                  double specified[{}], // specified = [{}]
                  double mass_matrix[{}], // computed
                  double forcing_vector[{}]); // computed
"""

    pyx_template = \
"""\
import numpy as np
cimport numpy as np

cdef extern from "{header_filename}":
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

    assert len(constants) == {mass_matrix_len}
    assert len(coordinates) == {coordinates_len}
    assert len(speeds) == {speeds_len}
    assert len(specified) == {specified_len}

    cdef np.ndarray[np.double_t, ndim=1, mode='c'] mass_matrix = np.zeros({mass_matrix_len})
    cdef np.ndarray[np.double_t, ndim=1, mode='c'] forcing_vector = np.zeros({forcing_vector_len})

    mass_forcing(<double*> constants.data,
                 <double*> coordinates.data,
                 <double*> speeds.data,
                 <double*> specified.data,
                 <double*> mass_matrix.data,
                 <double*> forcing_vector.data)

    return mass_matrix.reshape({forcing_vector_len}, {forcing_vector_len}), forcing_vector.reshape({forcing_vector_len}, 1)
"""
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

    # Comma separated lists of the input variables to the arrays, so that
    # you know which order you should supply them. These are used in the
    # comments in the C and header file.
    constant_list = ', '.join([str(c) for c in constants])
    coordinate_list = ', '.join([str(c).split('(')[0] for c in coordinates])
    speed_list = ', '.join([str(s).split('(')[0] for s in speeds])
    specified_list = ', '.join([str(e).split('(')[0] for e in specified])

    rows, cols = mass_matrix.shape
    mass_matrix_list = mass_matrix.reshape(rows * cols, 1).tolist()
    forcing_vector_list = forcing_vector.tolist()

    # TODO: make this optional
    sub_expressions, expressions = cse([entry[0] for entry in
                                        mass_matrix_list +
                                        forcing_vector_list],
                                       numbered_symbols('z_'))

    # make a dictionary that maps each symbol to the appropriate array index
    array_name_map = {'constants': constants,
                      'coordinates': coordinates,
                      'speeds': speeds,
                      'specified': specified}

    array_index_map = {}
    for array_name, variables in array_name_map.items():
        for i, var in enumerate(variables):
            array_index_map[str(var)] = r'{}[{}]'.format(array_name, i)

    def replace_with_array_index(string, replacement_map):
        unescaped_patterns = sorted(replacement_map.keys(), key=len,
                                    reverse=True)
        # escape parenatheses of variable names that are functions, i.e. have (t) in them
        escaped_patterns = [re.escape(k) for k in unescaped_patterns]
        pattern = re.compile('|'.join(escaped_patterns))

        def replacement(match):
            return replacement_map[match.group()]

        return pattern.sub(replacement, string)

    sub_expression_code_strings = []
    for var, exp in sub_expressions:
        code_str = replace_with_array_index(ccode(exp), array_index_map)
        sub_expression_code_strings.append('double {} = {};'.format(str(var), code_str))
    sub_expression_code_block = '\n    '.join(sub_expression_code_strings)

    mass_matrix_code_strings = []
    for i, exp in enumerate(expressions[:rows * cols]):
        code_str = replace_with_array_index(ccode(exp), array_index_map)
        mass_matrix_code_strings.append('{} = {};'.format('mass_matrix[{}]'.format(i), code_str))
    mass_matrix_code_block = '\n    '.join(mass_matrix_code_strings)

    forcing_vector_code_strings = []
    for i, exp in enumerate(expressions[rows * cols:]):
        code_str = replace_with_array_index(ccode(exp), array_index_map)
        forcing_vector_code_strings.append('{} = {};'.format('forcing_vector[{}]'.format(i), code_str))
    forcing_vector_code_block = '\n    '.join(forcing_vector_code_strings)

    c_code = c_template.format(header_filename, len(constants),
                               len(coordinates),
                               len(speeds), len(specified), rows
                               * cols, len(forcing_vector), constant_list,
                               coordinate_list, speed_list, specified_list,
                               sub_expression_code_block,
                               mass_matrix_code_block,
                               forcing_vector_code_block)

    with open(c_filename, 'w') as f:
        f.write(c_code)

    header_code = h_template.format(len(constants), constant_list,
                                    len(coordinates), coordinate_list,
                                    len(speeds), speed_list,
                                    len(specified), specified_list,
                                    rows * cols,
                                    len(forcing_vector))

    with open(header_filename, 'w') as f:
        f.write(header_code)

    pyx_code = pyx_template.format(header_filename=header_filename,
                                   mass_matrix_len=rows * cols,
                                   forcing_vector_len=len(forcing_vector),
                                   coordinates_len=len(coordinates),
                                   speeds_len=len(speeds),
                                   specified_len=len(speeds))

    with open(pyx_filename, 'w') as f:
        f.write(pyx_code)

    setup_code = setup_template.format(prefix=filename_prefix,
                                       c_filename=c_filename,
                                       pyx_filename=pyx_filename)

    with open(setup_py_filename, 'w') as f:
        f.write(setup_code)

def numeric_right_hand_side(kane, parameters, output_equations):
    """Returns the right hand side of the first order ordinary differential
    equations from a KanesMethod system which can be evaluated numerically.

    This function probably only works for simple cases (no dependent speeds,
    specified functions).

    Parameters
    ----------
    kane : KanesMethod object
    parameters : a sequence
        A list of symbols for the system constants.
    output_equations : list of expressions
        Each expression should be a function of the state and the parameters.

    Returns
    -------

    """

    # TODO: Deal with specified inputs.
    # TODO: Deal with problems with motion constraints.
    # TODO: Deal with problems with configuration constraints.
    # TODO: Should be able to handle the eom's in multiple forms: fr + frstar,
    # M & F, x' = f, etc

    dynamic = kane._q + kane._u
    dummy_symbols = [Dummy() for i in dynamic]
    dummy_dict = dict(zip(dynamic, dummy_symbols))
    kindiff_dict = kane.kindiffdict()

    M = kane.mass_matrix_full.subs(kindiff_dict).subs(dummy_dict)
    F = kane.forcing_full.subs(kindiff_dict).subs(dummy_dict)

    M_func = lambdify(dummy_symbols + parameters, M)
    F_func = lambdify(dummy_symbols + parameters, F)

    numerical_output_equations = []
    for equation in output_equations:
        dummified = equation.subs(kindiff_dict).subs(dummy_dict)
        numerical_output_equations.append(lambdify(dummy_symbols +
                                                   parameters, dummified))

    def right_hand_side(x, t, args):
        """Returns the derivatives of the states.

        Parameters
        ----------
        x : ndarray, shape(n)
            The current state vector.
        t : float
            The current time.
        args : ndarray
            The constants.

        Returns
        -------
        dx : ndarray, shape(2 * (n + 1))
            The derivative of the state.

        """
        arguments = np.hstack((x, args))
        dx = np.array(np.linalg.solve(M_func(*arguments),
            F_func(*arguments))).T[0]

        return dx

    def output(x, args):
        arguments = np.hstack((x, args))
        results = []
        for equation in numerical_output_equations:
            results.append(equation(arguments))

        return np.array(results)

    return right_hand_side, output
