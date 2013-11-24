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

# internal libraries
from templates import c_template, h_template, pyx_template, setup_template

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


def pydy_c_printer(constants, coordinates, speeds, specified=None):
    """Returns a subclass of sympy.printing.CCodPrinter to print appropriate
    C array index calls for all of the symbols in the equations of motion.

    For example:

    q(t) -> coordinates[1]
    m -> constants[2]
    F(t) -> specified[3]

    Parameters
    ----------
    constants : list of sympy.core.symbol.Symbol
        The constants in the equations of motion.
    coordinates : list of sympy.core.function.Function
        The generalized coordinates of the system.
    speeds : list of sympy.core.function.Function
        The generalized speeds of the system.
    specified : list of sympy.core.function.Function, optional, default=None
        The specifed quantities of the system.

    Returns
    -------
    PyDyCCodePrinter : class
        The modifed CCodePrinter class.

    """

    array_name_map = {'constants': constants,
                      'coordinates': coordinates,
                      'speeds': speeds}
    if specified is not None:
        array_name_map['specified'] = specified

    array_index_map = {}
    for array_name, variables in array_name_map.items():
        for i, var in enumerate(variables):
            array_index_map[var] = r'{}[{}]'.format(array_name, i)

    class PyDyCCodePrinter(CCodePrinter):

        def _print_Function(self, e):
            if e in array_index_map.keys():
                return array_index_map[e]
            else:
                return super(PyDyCCodePrinter, self)._print_Function(e)

        def _print_Symbol(self, e):
            if e in array_index_map.keys():
                return array_index_map[e]
            else:
                return super(PyDyCCodePrinter, self)._print_Symbol(e)

    return PyDyCCodePrinter


def generate_mass_forcing_cython_code(filename_prefix, mass_matrix,
                                      forcing_vector, constants,
                                      coordinates, speeds, specified=None):
    """Generates a Cython shared object module with a function that
    evaluates the mass_matrix and the forcing vector.

    Parameters
    ----------
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

    """

    # TODO : Maybe, allow expressions to be passed in for the specified
    # quantities for computation inside the c file. Would need the value of
    # time in the mass_forcing function.

    c_filename = filename_prefix + '_c.c'
    header_filename = filename_prefix + '_c.h'
    pyx_filename = filename_prefix + '.pyx'
    setup_py_filename = filename_prefix + '_setup.py'

    # Comma separated lists of the input variables to the arrays, so that
    # you know which order you should supply them. These are used in the
    # comments in the C and header file.
    # TODO: Add these to the doc string of the imported cythonized function.
    comma_lists = {}
    for list_name, sym_list in zip(['constants', 'coordinates', 'speeds',
                                    'specified'],
                                   [constants, coordinates, speeds,
                                    specified]):
        if sym_list is not None:
            comma_lists[list_name] = ', '.join([str(s).split('(')[0] for s
                                                in sym_list])
        else:
            comma_lists[list_name] = ['None Supplied']

    rows, cols = mass_matrix.shape
    mass_matrix_list = mass_matrix.reshape(rows * cols, 1).tolist()
    forcing_vector_list = forcing_vector.tolist()

    # TODO: make common sub expressions optional
    sub_expressions, expressions = cse([entry[0] for entry in
                                        mass_matrix_list +
                                        forcing_vector_list],
                                       numbered_symbols('z_'))

    expressions = {'common_sub': sub_expressions,
                   'mass_matrix': expressions[:rows * cols],
                   'forcing_vector': expressions[rows * cols:]}

    PyDyCCodePrinter = pydy_c_printer(constants, coordinates, speeds,
                                      specified)

    code_blocks = {}
    for exp_type, expression_list in expressions.items():
        c_lines = []
        for i, exp in enumerate(expression_list):
            if exp_type == 'common_sub':
                code_str = PyDyCCodePrinter().doprint(exp[1])
                c_lines.append('double {} = {};'.format(str(exp[0]),
                                                        code_str))
            else:
                code_str = PyDyCCodePrinter().doprint(exp)
                c_lines.append('{} = {};'.format('{}[{}]'.format(exp_type,
                                                                 i),
                                                 code_str))
        code_blocks[exp_type] = '\n    '.join(c_lines)

    template_values = {
        'header_filename': header_filename,
        'constants_len': len(constants),
        'mass_matrix_len': rows * cols,
        'forcing_vector_len': len(forcing_vector),
        'coordinates_len': len(coordinates),
        'speeds_len': len(speeds),
        'specified_len': len(specified) if specified is not None else 0,
        'constants_list': comma_lists['constants'],
        'coordinates_list': comma_lists['coordinates'],
        'speeds_list': comma_lists['speeds'],
        'specified_list': comma_lists['specified'],
        'sub_expression_block': code_blocks['common_sub'],
        'mass_matrix_block': code_blocks['mass_matrix'],
        'forcing_vector_block': code_blocks['forcing_vector'],
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
                                          specified=specified)

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
