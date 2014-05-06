#!/usr/bin/env python

# standard library
import os
import subprocess
import importlib
import random
from itertools import chain

# external libraries
import numpy as np
from sympy import lambdify, numbered_symbols, cse, symbols
from sympy.printing.ccode import CCodePrinter
try:
    import theano
except ImportError:  # If theano is not installed.
    theano_installed = False
else:
    from sympy.printing.theanocode import theano_function
    theano_installed = True
try:
    import Cython
except ImportError:  # If Cython is not installed.
    cython_installed = False
else:
    cython_installed = True

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


class CythonGenerator(object):

    def __init__(self, filename_prefix, mass_matrix, forcing_vector,
                 constants, coordinates, speeds, specified=None):

        """Instantiates an object that can generates a Cython shared object
        module with a function that evaluates the provided mass_matrix and
        the forcing vector given the numerical values of the input
        variables.

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

        self.filename_prefix = filename_prefix
        self.mass_matrix = mass_matrix
        self.forcing_vector = forcing_vector
        self.symbols = {'constants': constants,
                        'coordinates': coordinates,
                        'speeds': speeds,
                        'specified': specified}

        self.rows, self.cols = mass_matrix.shape
        self.expressions = \
            {'mass_matrix': mass_matrix.reshape(self.rows * self.cols,
                                                1).tolist(),
             'forcing_vector': forcing_vector.tolist()}

        self._generate_file_names()
        self._generate_comma_lists()
        self._generate_sub_expressions()
        self._generate_pydy_c_printer()
        self._generate_code_blocks()
        self._create_template_dict()

    def _generate_file_names(self):
        """Generates the names for the four files needed to build the Cython
        extension module."""

        self.c_filename = self.filename_prefix + '_c.c'
        self.header_filename = self.filename_prefix + '_c.h'
        self.pyx_filename = self.filename_prefix + '.pyx'
        self.setup_py_filename = self.filename_prefix + '_setup.py'

    def _generate_comma_lists(self):
        """Generates comma separated lists of the input variables to the
        arrays, so that you know which order you should supply them. These
        are used in the comments in the C and header file."""

        # TODO: Add these to the doc string of the imported cythonized
        # function.

        self.comma_lists = {}
        for list_name, sym_list in self.symbols.items():
            if sym_list is not None:
                self.comma_lists[list_name] = \
                    ', '.join([str(s).split('(')[0] for s in sym_list])
            else:
                self.comma_lists[list_name] = ['None Supplied']

    def _generate_sub_expressions(self, int_var='z_'):
        """Finds common subexpressions in the mass matrix and forcing vector
        expression lists and rewrites them."""

        # TODO: make common sub expressions optional

        list_of_lists = self.symbols.values()
        if None in list_of_lists:
            list_of_lists.remove(None)
        all_symbols = list(chain.from_iterable(list_of_lists))
        while symbols(int_var) in all_symbols:
            int_var = random.choice(all_letters) + '_'

        sub_expressions, expressions = \
            cse([entry[0] for entry in self.expressions['mass_matrix'] +
                 self.expressions['forcing_vector']],
                numbered_symbols(int_var))

        self.expressions = \
            {'common_sub': sub_expressions,
             'mass_matrix': expressions[:self.rows * self.cols],
             'forcing_vector': expressions[self.rows * self.cols:]}

    def _generate_pydy_c_printer(self):
        """Returns a subclass of sympy.printing.CCodePrinter to print
        appropriate C array index calls for all of the symbols in the equations
        of motion.

        Examples
        --------

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> from pydy.codegen.code import CythonGenerator
        >>> cython_generator = CythonGenerator(...)
        >>> cython_generator._generate_pydy_c_printer()
        >>> m = symbols('m') # m is the first constant in the EoMs
        >>> cython_generator.PyDyCCodePrinter().doprint(m)
        constants[0];
        >>> q = dynamicsymbols('q') # q is the second coordinate in the EoMs
        >>> cython_generator.PyDyCCodePrinter().doprint(q)
        coordinates[1];
        >>> F = dynamicsymbols('F') # F is the third specified in the EoMs
        >>> cython_generator.PyDyCCodePrinter().doprint(F)
        specified[2];

        """

        array_index_map = {}
        for array_name, variables in self.symbols.items():
            if variables is not None:
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

        self.PyDyCCodePrinter = PyDyCCodePrinter

    def _generate_code_blocks(self):
        """Writes the blocks of code for the C file."""

        self.code_blocks = {}
        for exp_type, expression_list in self.expressions.items():
            c_lines = []
            for i, exp in enumerate(expression_list):
                if exp_type == 'common_sub':
                    code_str = self.PyDyCCodePrinter().doprint(exp[1])
                    c_lines.append('double {} = {};'.format(str(exp[0]),
                                                            code_str))
                else:
                    code_str = self.PyDyCCodePrinter().doprint(exp)
                    lhs = '{}[{}]'.format(exp_type, i)
                    c_lines.append('{} = {};'.format(lhs, code_str))
            self.code_blocks[exp_type] = '\n    '.join(c_lines)

    def _create_template_dict(self):
        """Creates a dictionary mapping all of the variables in the
        templates to the appropriate strings."""

        self.template_values = {
            'header_filename': self.header_filename,
            'mass_matrix_len': self.rows * self.cols,
            'forcing_vector_len': len(self.forcing_vector),
            'constants_len': len(self.symbols['constants']),
            'coordinates_len': len(self.symbols['coordinates']),
            'speeds_len': len(self.symbols['speeds']),
            'specified_len': (len(self.symbols['specified']) if
                              self.symbols['specified'] is not None else 0),
            'constants_list': self.comma_lists['constants'],
            'coordinates_list': self.comma_lists['coordinates'],
            'speeds_list': self.comma_lists['speeds'],
            'specified_list': self.comma_lists['specified'],
            'sub_expression_block': self.code_blocks['common_sub'],
            'mass_matrix_block': self.code_blocks['mass_matrix'],
            'forcing_vector_block': self.code_blocks['forcing_vector'],
            'prefix': self.filename_prefix,
            'c_filename': self.c_filename,
            'pyx_filename': self.pyx_filename,
        }

        if self.symbols['specified'] is not None:
            specified_template = {
                'specified_double': " " * 18 + "double specified[{specified_len}], // specified = [{specified_list}]".format(**self.template_values) + '\n',
                'specified_assert': "\n    assert len(specified) == {specified_len}\n".format(**self.template_values),
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

        self.template_values.update(specified_template)

    def _write_cython_code(self):
        """Writes the files needed for the Cython extensions module in the
        current directory."""

        files = {self.c_filename: c_template,
                 self.header_filename: h_template,
                 self.pyx_filename: pyx_template,
                 self.setup_py_filename: setup_template}

        for filename, template in files.items():
            code = template.format(**self.template_values)
            with open(filename, 'w') as f:
                f.write(code)

    def _compile_cython_code(self):
        """Compiles the Cython extension module using distutils."""

        if not cython_installed:
            raise ValueError('Cython is not installed.')

        # TODO : Need some way to cleanup the files creates by this after
        # use.

        # TODO : This may not be cross platform. Needs to be explored on
        # Windows and Mac.

        # This prevents output to stdout and waits till it is done.
        cmd = ['python', self.setup_py_filename, 'build_ext', '--inplace']
        subprocess.call(cmd, stderr=subprocess.STDOUT,
                        stdout=subprocess.PIPE)

    def generate_extension(self):
        """Generates a Cython extensions module with the given file name
        prefix which contains a function `mass_forcing_matrices` that
        evaluates the mass matrix and forcing function."""
        self._write_cython_code()
        self._compile_cython_code()


def generate_ode_function(mass_matrix, forcing_vector, constants,
                          coordinates, speeds, specified=None,
                          generator='lambdify'):
    """Returns a numerical function which can evaluate the right hand side
    of the first order ordinary differential equations from a system
    described by:

    M(constants, coordinates) x' = F(constants, coordinates, speeds, specified)

    Parameters
    ----------
    mass_matrix : sympy.Matrix, shape(n,n)
        The symbolic mass matrix of the system.
    forcing_vector : sympy.Matrix, shape(n,1)
        The symbolic forcing vector of the system.
    constants : list of sympy.Symbol
        The constants in the equations of motion.
    coordinates : list of sympy.Function
        The generalized coordinates of the system.
    speeds : list of sympy.Function
        The generalized speeds of the system.
    specified : list of sympy.Function
        The specifed quantities of the system.
    generator : string, {'lambdify'|'theano'|'cython'}, optional
        The method used for generating the numeric right hand side.

    Returns
    -------
    evaluate_ode_function : function
        A function which evaluates the derivaties of the states.

    """
    if generator == 'theano' and not theano_installed:
        raise ValueError('Theano is not installed.')

    if generator == 'cython' and not cython_installed:
        raise ValueError('Cython is not installed.')

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
        else:
            raise ImportError('Theano is not installed, choose another method.')

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

            return (mass_matrix_func(*value_array),
                    forcing_vector_func(*value_array))

    elif generator == 'cython':

        filename_prefix = 'multibody_system'

        # TODO : This is a hack to allow you to regenerate cython modules
        # without closing the Python session. It may be best to also force
        # the user to provide a module name when generating the Cython code.
        # Check out the Cython inline code to figure out how to do all this
        # better with disutils:
        # https://github.com/cython/cython/blob/master/Cython/Build/Inline.py

        # The .pyx file has the same prefix as the Cython generated [.dll,
        # .so, .dylib] shared library file, so we should be able to check
        # all files in the directory for matches except the .pyx file.
        prefixes = [os.path.splitext(p)[0] for p in os.listdir('.') if not
                    p.endswith('.pyx')]
        while True:
            if filename_prefix in prefixes:
                filename_prefix += '_' + random.choice(all_letters)
            else:
                break

        cython_generator = CythonGenerator(filename_prefix, mass_matrix,
                                           forcing_vector, constants,
                                           coordinates, speeds,
                                           specified=specified)
        cython_generator.generate_extension()

        cython_module = importlib.import_module(filename_prefix)
        mass_forcing_func = cython_module.mass_forcing_matrices

    else:
        # TODO : add numba, fortran, parakeet, sympy.autowrap (needs matrix
        # support)
        raise NotImplementedError('The {} code generation is not implemented'.format(generator))

    def evaluate_ode(x, t, args):
        """Returns the derivatives of the states, i.e. numerically evaluates
        the right hand side of the first order differential equation(s).

        x' = f(x, t)

        Parameters
        ----------
        x : ndarray, shape({num_states},)
            The current state vector:
                {state_list}
        t : float
            The current time.
        args : dictionary
            constants : ndarray, shape({num_constants},)
                {constant_list}
            specified : ndarray, shape({num_specified},) or a function
                If this is a function it must be of the form f(x, t), where
                x is the current state vector and t is the current time and
                it must return an ndarray of the correct shape.
                {specified_list}

        Returns
        -------
        dx : ndarray, shape({num_states},)
            The derivative of the state vector.

        """

        segmented = [args['constants'],
                     x[:len(coordinates)],
                     x[len(coordinates):]]

        if specified is not None:
            try:
                sp_val = args['specified'](x, t)
            except TypeError:  # not callable
                # If not callable, then it should be a float or ndarray.
                sp_val = args['specified']

            # If the value is just a float, then convert to a 1D array.
            try:
                len(sp_val)
            except TypeError:
                sp_val = np.asarray([sp_val])

            segmented.append(sp_val)

        mass_matrix_values, forcing_vector_values = \
            mass_forcing_func(*segmented)

        # TODO: figure out how to off load solve to the various generated
        # code, for example for Theano:
        # http://deeplearning.net/software/theano/library/sandbox/linalg.html#theano.sandbox.linalg.ops.Solve

        # Could use scipy.linalg.solve and enable a and b overwriting to
        # avoid the array copying.
        dx = np.array(np.linalg.solve(mass_matrix_values,
                                      forcing_vector_values)).T[0]

        return dx

    template_values = {'num_states': len(coordinates + speeds),
                       'state_list': ', '.join([str(s) for s in coordinates
                                                + speeds]),
                       'num_constants': len(constants),
                       'constant_list': ', '.join([str(c) for c in constants]),
                       'num_specified': '0',
                       'specified_list': '',
                       }

    if specified is not None:
        template_values['num_specified'] = len(specified)
        template_values['specified_list'] = ', '.join([str(s) for s in
                                                       specified])

    evaluate_ode.__doc__ = evaluate_ode.__doc__.format(**template_values)

    return evaluate_ode
