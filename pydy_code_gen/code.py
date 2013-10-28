#!/usr/bin/env python

# standard library
import time
import re

# external libraries
import numpy as np
from sympy import Dummy, lambdify, numbered_symbols, cse
from sympy.printing import ccode
from sympy.utilities.autowrap import autowrap
from sympy.printing.theanocode import theano_function

# debugging
try:
    from IPython.core.debugger import Tracer
except ImportError:
    pass
else:
    set_trace = Tracer()


def generate_numeric_eom_matrices(mass_matrix, forcing_vector, constants,
                                  states, specified=None,
                                  generator="lambdify", **options):
    """Returns functions which compute the mass matrix and forcing vector
    given numerical values of the states, constants, and optionally
    specified variables.

    Parameters
    ----------
    mass_matrix : SymPy Matrix, shape(n, n)
        An n x n symbolic matrix where each entry is an expression made up of
        the states, constants, and specified variables.
    forcing_vector : SymPy Matrix, shape(n,)
        An n x 1 symbolic matrix where each entry is an expression made up of
        the states, constants, and specified variables.
    constants : sequence of sympy.Symbol
    states : sequence of sympy.Function
    specified : sequence of sympy.Function
    generator : string, optional, default='lambdify'

    Returns
    -------
    mass_matrix_func : function
    forcing_vector_func : function

    constants + states [+ specified]

    Examples
    --------

    >>> q = dynamicsymbols('q:2')
    >>> u = dynamicsymbols('u:2')
    >>> f = dynamicsymbols('f:2')
    >>> c = symbols('c:2')
    >>> mass_matrix = Matrix([[q[0] + u[0] + c[0], c[0] * (q[1] + u[1])],

    >>> type(mass_matrix)
    sympy.matrices.dense.MutableDenseMatrix
    >>> type(forcing_vector)
    sympy.matrices.dense.MutableDenseMatrix
    >>> generate_numeric_matrices(mass_matrix, forcing_vector, c, q + u,
        >>> specified=f,


    """

    dynamic = states
    if specified is not None:
        dynamic += specified

    # TODO : the Dummy symbols may not be needed now that lambdify is updated.
    # the theano option may not need them either.
    dummy_symbols = [Dummy() for i in dynamic]
    dummy_dict = dict(zip(dynamic, dummy_symbols))

    dummy_mass_matrix = mass_matrix.subs(dummy_dict)
    dummy_forcing_vector = forcing_vector.subs(dummy_dict)

    arguments = constants + dummy_symbols

    if generator == 'lambdify':
        mass_matrix_func = lambdify(arguments, dummy_mass_matrix)
        forcing_vector_func = lambdify(arguments, dummy_forcing_vector)
    elif generator == 'theano':
        mass_matrix_func = theano_function(arguments, [dummy_mass_matrix],
                                           on_unused_input='ignore')
        forcing_vector_func = theano_function(arguments,
                                              [dummy_forcing_vector],
                                              on_unused_input='ignore')
        # lower run time from 0.15s to 0.08s for n=1
        mass_matrix_func.trust_input = True
        forcing_vector_func.trust_input = True

    elif generator == 'autowrap':
        funcs = []
        for entry in dummy_mass_matrix:
            funcs.append(autowrap(entry, args=arguments, **options))

        def mass_matrix_func(*args):
            result = []
            for func in funcs:
                result.append(func(*args))
            # TODO : this may not be correctly reshaped
            return np.matrix(result).reshape(np.sqrt(len(result)),
                                             np.sqrt(len(result)))

        funcs = []
        for row in dummy_forcing_vector:
            funcs.append(autowrap(row, args=arguments, **options))

        def forcing_vector_func(*args):
            result = []
            for func in funcs:
                result.append(func(*args))
            return np.matrix(result)
    else:
        raise NotImplementedError('{} is not implemented yet'.format(generator))

    return mass_matrix_func, forcing_vector_func


def generate_mass_forcing_cython_code(filename_prefix, mass_matrix,
                                      forcing_vector, constants,
                                      coordinates, speeds, specified=None,
                                      time_variable='t'):

    # TODO : Maybe, allow expressions to be passed in for the specified
    # quantities for computation inside the c file. We would need the value
    # of time in the mass_forcing function.

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
    # TODO: Add this to the doc string of the imported cythonized function.
    constant_list = ', '.join([str(c) for c in constants])
    coordinate_list = ', '.join([str(c).split('(')[0] for c in coordinates])
    speed_list = ', '.join([str(s).split('(')[0] for s in speeds])
    if specified is not None:
        specified_list = ', '.join([str(e).split('(')[0] for e in specified])
    else:
        specified_list = ['None Supplied'] # This will cause some issues if the person has variable named the same.

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
        # escape parentheses of variable names that are functions, i.e. have (t) in them
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

    c_code = c_template.format(header_filename,
                               len(constants),
                               len(coordinates),
                               len(speeds),
                               len(specified),
                               rows * cols,
                               len(forcing_vector),
                               constant_list,
                               coordinate_list,
                               speed_list,
                               specified_list,
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
                                   specified_len=len(specified))

    with open(pyx_filename, 'w') as f:
        f.write(pyx_code)

    setup_code = setup_template.format(prefix=filename_prefix,
                                       c_filename=c_filename,
                                       pyx_filename=pyx_filename)

    with open(setup_py_filename, 'w') as f:
        f.write(setup_code)

def numeric_right_hand_side(kane, parameters, specified=None, generator='lambdify'):
    """Returns the right hand side of the first order ordinary differential
    equations from a KanesMethod system which can be evaluated numerically.

    This function probably only works for simple cases (no dependent speeds,
    specified functions, etc).

    Parameters
    ----------
    kane : KanesMethod object
    parameters : a sequence
        A list of symbols for the system constants.
    specified : a sequence
        A sequence of symbols/functions for specifed variables.
    generator : optional
        lambdify theano numba cython fortran parakeet

    Returns
    -------

    """

    dynamic = kane._q + kane._u

    kindiff_dict = kane.kindiffdict()

    mass_matrix = kane.mass_matrix_full.subs(kindiff_dict)
    forcing_vector = kane.forcing_full.subs(kindiff_dict)

    mass_matrix_func, forcing_vector_func = \
        generate_numeric_eom_matrices(mass_matrix, forcing_vector,
                                      parameters, dynamic,
                                      specified=specified,
                                      generator=generator)

    arguments = parameters + dynamic
    if specified is not None: arguments += specified
    start = time.time()
    numtimes = 1000

    # Lower from 0.5s to 0.15s for the run time when n=1
    #inp = np.random.random(len(arguments))
    inp = [np.asarray(x) for x in np.random.random(len(arguments))]
    for i in range(numtimes):
        mass_matrix_func(*inp)
        forcing_vector_func(*inp)
    total = time.time() - start
    print("It took {} seconds to compute M and F with {} {} times at an \
average of {} seconds per computation.".format(total, generator, numtimes,
    total / numtimes))

    #set_trace()

    def right_hand_side(x, t, args):
        """Returns the derivatives of the states.

        Parameters
        ----------
        x : ndarray, shape(n)
            The current state vector.
        t : float
            The current time.
        args : ndarray
            The specified variables and constants.

        Returns
        -------
        dx : ndarray, shape(2 * (n + 1))
            The derivative of the state.

        """
        arguments = [np.asarray(x) for x in np.hstack((x, args))]
        # TODO: figure out how to off load solve to Theano
        # http://deeplearning.net/software/theano/library/sandbox/linalg.html#theano.sandbox.linalg.ops.Solve
        dx = np.array(np.linalg.solve(mass_matrix_func(*arguments),
                                      forcing_vector_func(*arguments))).T[0]

        return dx

    return right_hand_side
