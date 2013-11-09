#!/usr/bin/env python

# standard library
import subprocess
import importlib
import random

# external libraries
import numpy as np
from sympy import lambdify, numbered_symbols, cse
from sympy.printing.ccode import CCodePrinter
from sympy.printing.theanocode import theano_function

# Python 2 vs 3 importing
try:
    from string import letters as all_letters
except ImportError:
    from string import ascii_letters as all_letters

# debugging
try:
    from IPython.core.debugger import Tracer
except ImportError:
    pass
else:
    set_trace = Tracer()


def generate_mass_forcing_cython_code(filename_prefix, mass_matrix,
                                      forcing_vector, constants,
                                      coordinates, speeds, specified=None,
                                      time_variable='t'):
    """Generates a Cython shared object module with a function that
    evaluates the mass_matrix and the forcing vector.

    Parameters
    ==========
    filename_prefix : string
        The desired name of the created module.
    mass_matrix : sympy.matrices.dense.MutableDenseMatrix, shape(n,n)
        The symbolic mass matrix of the system.
    forcing_vector : sympy.matrices.dense.MutableDenseMatrix, shape(n,1)
        The symbolic forcing vector of the system.
    constants : list of sympy.core.symbol.Symbol
        The constants in the equations of motion.
    coordinates : list of sympy.core.function.Function
        The generalized coordinates of the system.
    speeds : list of sympy.core.function.Function
        The generalized speeds of the system.
    specified : list of sympy.core.function.Function, optional, default=None
        The specifed quantities of the system.
    time_variable : str, option, default='t'
        If you use a varialbe name for time other than the default in
        sympy.physics.mechanics, you will need to specifiy that here.

    """

    # TODO : Maybe, allow expressions to be passed in for the specified
    # quantities for computation inside the c file. Would need the value of
    # time in the mass_forcing function.

    c_filename = filename_prefix + '_c.c'
    header_filename = filename_prefix + '_c.h'
    pyx_filename = filename_prefix + '.pyx'
    setup_py_filename = filename_prefix + '_setup.py'

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
}}"""

    h_template = \
"""\
void mass_forcing(double constants[{constants_len}], // constants = [{constants_list}]
                  double coordinates[{coordinates_len}], // coordinates = [{coordinates_list}]
                  double speeds[{speeds_len}], // speeds = [{speeds_list}]
{specified_double}                  double mass_matrix[{mass_matrix_len}], // computed
                  double forcing_vector[{forcing_vector_len}]); // computed"""

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
                      'speeds': speeds}
    if specified is not None:
        array_name_map['specified'] = specified

    array_index_map = {}
    for array_name, variables in array_name_map.items():
        for i, var in enumerate(variables):
            array_index_map[str(var)] = r'{}[{}]'.format(array_name, i)

    class pydy_ccode(CCodePrinter):
        def _print_Function(self, e):
            if str(e) in array_index_map.keys():
                return array_index_map[str(e)]
            else:
                return super()._print_Function(e)
        def _print_Symbol(self, e):
            if str(e) in array_index_map.keys():
                return array_index_map[str(e)]
            else:
                return super()._print_Symbol(e)

    sub_expression_code_strings = []
    for var, exp in sub_expressions:
        code_str = pydy_ccode().doprint(exp)
        sub_expression_code_strings.append('double {} = {};'.format(str(var), code_str))
    sub_expression_code_block = '\n    '.join(sub_expression_code_strings)

    mass_matrix_code_strings = []
    for i, exp in enumerate(expressions[:rows * cols]):
        code_str = pydy_ccode().doprint(exp)
        mass_matrix_code_strings.append('{} = {};'.format('mass_matrix[{}]'.format(i), code_str))
    mass_matrix_code_block = '\n    '.join(mass_matrix_code_strings)

    forcing_vector_code_strings = []
    for i, exp in enumerate(expressions[rows * cols:]):
        code_str = pydy_ccode().doprint(exp)
        forcing_vector_code_strings.append('{} = {};'.format('forcing_vector[{}]'.format(i), code_str))
    forcing_vector_code_block = '\n    '.join(forcing_vector_code_strings)

    template_values = {
        'header_filename': header_filename,
        'constants_len': len(constants),
        'mass_matrix_len': rows * cols,
        'forcing_vector_len': len(forcing_vector),
        'coordinates_len': len(coordinates),
        'speeds_len': len(speeds),
        'specified_len': len(specified) if specified is not None else 0,
        'constants_list': constant_list,
        'coordinates_list': coordinate_list,
        'speeds_list': speed_list,
        'specified_list': specified_list,
        'sub_expression_block': sub_expression_code_block,
        'mass_matrix_block': mass_matrix_code_block,
        'forcing_vector_block': forcing_vector_code_block,
        'prefix': filename_prefix,
        'c_filename': c_filename,
        'pyx_filename': pyx_filename,
    }

    if specified is not None:
        specified_template = {
            'specified_double': " " * 18 + "double specified[{specified_len}], // specified = [{specified_list}]".format(**template_values) + '\n',
            'specified_assert': "\n    assert len(specified) == {specified_len}\n".format(**template_values),
            'def_specified_arg': ",\n                          np.ndarray[np.double_t, ndim=1, mode='c'] specified",
            'cdef_specified_arg': "\n" + ' ' * 22 + "double* specified,",
            'call_specified_arg': "\n" + ' ' * 17 + "<double*> specified.data,"
        }
    else:
        specified_template = {
            'specified_double': "",
            'specified_assert': "",
            'def_specified_arg': "",
            'cdef_specified_arg': "",
            'call_specified_arg': ""
        }

    template_values.update(specified_template)

    files = {c_filename: c_template,
             header_filename: h_template,
             pyx_filename: pyx_template,
             setup_py_filename: setup_template}

    for filename, template in files.items():

        code = template.format(**template_values)

        with open(filename, 'w') as f:
            f.write(code)

    # This prevents output to stdout and waits till it is done.
    p = subprocess.Popen(['python', setup_py_filename, 'build_ext',
                          '--inplace'], stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)
    p.wait()

    # TODO : Need some way to cleanup the files creates by this after use.


def numeric_right_hand_side(mass_matrix, forcing_vector, constants,
                            coordinates, speeds, specified=None,
                            generator='lambdify'):
    """Returns a function for the right hand side of the first order
    ordinary differential equations from a system described by:

    M(constants, coordinates) x' = F(constants, coordinates, speeds, specified)

    which can be evaluated numerically.

    Parameters
    ----------
    mass_matrix : sympy.matrices.dense.MutableDenseMatrix, shape(n,n)
        The symbolic mass matrix of the system.
    forcing_vector : sympy.matrices.dense.MutableDenseMatrix, shape(n,1)
        The symbolic forcing vector of the system.
    constants : list of sympy.core.symbol.Symbol
        The constants in the equations of motion.
    coordinates : list of sympy.core.function.Function
        The generalized coordinates of the system.
    speeds : list of sympy.core.function.Function
        The generalized speeds of the system.
    specified : list of sympy.core.function.Function
        The specifed quantities of the system.
    generator : string, {'lambdify'|'theano'|'cython'}, optional
        The method used for generating the numeric right hand side.

    Returns
    -------
    right_hand_side : function
        A function which evaluates the derivaties of the states.

    """

    if generator == 'lambdify' or generator == 'theano':

        arguments = constants + coordinates + speeds
        if specified is not None:
            arguments += specified

        if generator == 'lambdify':

            mass_matrix_func = lambdify(arguments, mass_matrix)
            forcing_vector_func = lambdify(arguments, forcing_vector)

        elif generator == 'theano':

            mass_matrix_func = theano_function(arguments, [mass_matrix],
                                               on_unused_input='ignore')
            forcing_vector_func = theano_function(arguments,
                                                  [forcing_vector],
                                                  on_unused_input='ignore')
            # Theano will run faster if you trust the input. I'm not sure
            # what the implications of this are. See:
            # http://deeplearning.net/software/theano/tutorial/faq.html#faster-small-theano-function
            mass_matrix_func.trust_input = True
            forcing_vector_func.trust_input = True

        def mass_forcing_func(numerical_constants, numerical_coordinates,
                              numerical_speeds, numerical_specified=None):
            """Returns numerical evaluations of the mass matrix and forcing
            vector."""

            values = [numerical_constants, numerical_coordinates,
                      numerical_speeds]
            if specified is not None:
                values.append(numerical_specified)

            value_array = np.hstack(tuple(values))
            if generator == 'theano':
                value_array = [np.asarray(v) for v in value_array]

            return mass_matrix_func(*value_array), forcing_vector_func(*value_array)

    elif generator == 'cython':

        filename_prefix = 'multibody_system'

        # TODO : This is a hack to allow you to regenerate cython modules
        # without closing the Python session. It may be best to also force
        # the user to provide a module name when generating the Cython code.
        # Check out the Cython inline code to figure out how to do all this
        # better with disutils:
        # https://github.com/cython/cython/blob/master/Cython/Build/Inline.py
        exists = True
        while exists:
            try:
                open(filename_prefix + '.so', 'r')
            except IOError:
                exists = False
            else:
                filename_prefix += '_' + random.choice(all_letters)

        generate_mass_forcing_cython_code(filename_prefix, mass_matrix,
                                          forcing_vector, constants,
                                          coordinates, speeds,
                                          specified=specified,
                                          time_variable='t')

        cython_module = importlib.import_module(filename_prefix)
        mass_forcing_func = cython_module.mass_forcing_matrices

    else:
        # TODO : add numba, fortan, parakeet, sympy.autowrap (needs matrix
        # support)
        raise NotImplementedError('The {} code generation is not implemented'.format(generator))

    def right_hand_side(x, t, args):
        """Returns the derivatives of the states.

        Parameters
        ----------
        x : ndarray, shape(n,)
            The current state vector.
        t : float
            The current time.
        args : dictionary
            constants : ndarray
            specified : ndarray
            num_coordinates : integer

        Returns
        -------
        dx : ndarray, shape(n,)
            The derivative of the state.

        """
        # TODO : Add more useful info to this doc string. Generate it
        # dynamically.
        # http://stackoverflow.com/questions/10307696/how-to-put-a-variable-into-python-docstring

        # TODO : Allow arg['specified'] to be a function of the states and
        # time, which gets evaluated at each time step.

        segmented = [args['constants'], x[:args['num_coordinates']],
                     x[args['num_coordinates']:]]
        if specified is not None:
            segmented.append(args['specified'])

        mass_matrix_values, forcing_vector_values = \
            mass_forcing_func(*segmented)

        # TODO: figure out how to off load solve to the various generated
        # code, for example for Theano:
        # http://deeplearning.net/software/theano/library/sandbox/linalg.html#theano.sandbox.linalg.ops.Solve

        dx = np.array(np.linalg.solve(mass_matrix_values,
                                      forcing_vector_values)).T[0]

        return dx

    return right_hand_side
