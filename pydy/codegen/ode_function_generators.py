#!/usr/bin/env python

import sys
if sys.version_info > (3, 0):
    from collections.abc import Sequence
else:
    from collections import Sequence
from itertools import chain

import numpy as np
import numpy.linalg
import scipy.linalg
import sympy as sm
import sympy.physics.mechanics as me
from sympy.core.function import UndefinedFunction, Derivative
Cython = sm.external.import_module('Cython')
theano = sm.external.import_module('theano')
if theano:
    from sympy.printing.theanocode import theano_function

from .c_code import _CLUsolveGenerator
from .cython_code import CythonMatrixGenerator


class ODEFunctionGenerator(object):
    """This is an abstract base class for all of the generators. A subclass
    is expected to implement the methods necessary to evaluate the arrays
    needed to compute xdot for the three different system specification
    types."""

    _rhs_doc_template = \
"""\
Returns the derivatives of the states, i.e. numerically evaluates the right
hand side of the first order differential equation.

x' = f(x, t,{specified_call_sig} p)

Parameters
==========
x : ndarray, shape({num_states},)
    The state vector is ordered as such:
{state_list}
t : float
    The current time.{specifieds_explanation}{constants_explanation}

Returns
=======
dx : ndarray, shape({num_states},)
    The derivative of the state vector.

"""

    _constants_doc_templates = {}

    _constants_doc_templates[None] = \
"""
p : dictionary len({num_constants}) or ndarray shape({num_constants},)
    Either a dictionary that maps the constants symbols to their numerical
    values or an array with the constants in the following order:
{constant_list}\
"""

    _constants_doc_templates['array'] = \
"""
p : ndarray shape({num_constants},)
    A ndarray of floats that give the numerical values of the constants in
    this order:
{constant_list}\
"""

    _constants_doc_templates['dictionary'] = \
"""
p : dictionary len({num_constants})
    A dictionary that maps the constants symbols to their numerical values
    with at least these keys:
{constant_list}\
"""

    _specifieds_doc_templates = {}

    _specifieds_doc_templates[None] = \
"""
r : dictionary; ndarray, shape({num_specified},); function

    There are three options for this argument. (1) is more flexible but
    (2) and (3) are much more efficient.

    (1) A dictionary that maps the specified functions of time to floats,
    ndarrays, or functions that produce ndarrays. The keys can be a single
    specified symbolic function of time or a tuple of symbols. The total
    number of symbols must be equal to {num_specified}. If the value is a
    function it must be of the form g(x, t), where x is the current state
    vector ndarray and t is the current time float and it must return an
    ndarray of the correct shape. For example::

      r = {{a: 1.0,
           (d, b) : np.array([1.0, 2.0]),
           (e, f) : lambda x, t: np.array(x[0], x[1]),
           c: lambda x, t: np.array(x[2])}}

    (2) A ndarray with the specified values in the correct order and of the
    correct shape.

    (3) A function that must be of the form g(x, t), where x is the current
    state vector and t is the current time and it must return an ndarray of
    the correct shape.

    The specified inputs are, in order:
{specified_list}\
"""

    _specifieds_doc_templates['array'] = \
"""
r : ndarray, shape({num_specified},)

    A ndarray with the specified values in the correct order and of the
    correct shape.

    The specified inputs are, in order:
{specified_list}\
"""

    _specifieds_doc_templates['function'] = \
"""
r : function

    A function that must be of the form g(x, t), where x is the current
    state vector and t is the current time and it must return an ndarray of
    shape({num_specified},).

    The specified inputs are, in order:
{specified_list}\
"""

    _specifieds_doc_templates['dictionary'] = \
"""
r : dictionary

    A dictionary that maps the specified functions of time to floats,
    ndarrays, or functions that produce ndarrays. The keys can be a single
    specified symbolic function of time or a tuple of symbols. The total
    number of symbols must be equal to {num_specified}. If the value is a
    function it must be of the form g(x, t), where x is the current state
    vector ndarray and t is the current time float and it must return an
    ndarray of the correct shape. For example::

      r = {{a: 1.0,
           (d, b) : np.array([1.0, 2.0]),
           (e, f) : lambda x, t: np.array(x[0], x[1]),
           c: lambda x, t: np.array(x[2])}}

    The specified inputs are, in order:
{specified_list}\
"""

    @staticmethod
    def _deduce_system_type(**kwargs):
        """Based on the combination of arguments this returns which ODE
        description has been provided.

        full rhs
            x' = f(x, t, r, p)
        full mass matrix
            M(x, p) * x' = f(x, t, r, p)
        min mass matrix
            M(q, p) * u' = f(q, u, t, r, p)
            q' = g(q, u, t)

        """

        if kwargs.pop('coordinate_derivatives') is not None:
            system_type = 'min mass matrix'
        elif kwargs.pop('mass_matrix') is not None:
            system_type = 'full mass matrix'
        else:
            system_type = 'full rhs'

        return system_type

    def __init__(self, right_hand_side, coordinates, speeds, constants=(),
                 mass_matrix=None, coordinate_derivatives=None,
                 specifieds=None, linear_sys_solver='numpy',
                 constants_arg_type=None, specifieds_arg_type=None):
        """Generates a numerical function which can evaluate the right hand
        side of the first order ordinary differential equations from a
        system described by one of the following three symbolic forms:

            [1] x' = F(x, t, r, p)

            [2] M(x, p) x' = F(x, t, r, p)

            [3] M(q, p) u' = F(q, u, t, r, p)
                q' = G(q, u, t, r, p)

        where

            x : states, i.e. [q, u]
            t : time
            r : specified (exogenous) inputs
            p : constants
            q : generalized coordinates
            u : generalized speeds
            M : mass matrix (full or minimum)
            F : right hand side (full or minimum)
            G : right hand side of the kinematical differential equations

        The generated function is of the form F(x, t, p) or F(x, t, r, p)
        depending on whether the system has specified inputs or not.

        Parameters
        ==========
        right_hand_side : SymPy Matrix, shape(n, 1)
            A column vector containing the symbolic expressions for the
            right hand side of the ordinary differential equations. If the
            right hand side has been solved for symbolically then only F is
            required, see form [1]; if not then the mass matrix must also be
            supplied, see forms [2, 3].
        coordinates : sequence of SymPy Functions
            The generalized coordinates. These must be ordered in the same
            order as the rows in M, F, and/or G and be functions of time.
        speeds : sequence of SymPy Functions
            The generalized speeds. These must be ordered in the same order
            as the rows in M, F, and/or G and be functions of time.
        constants : sequence of SymPy Symbols, optional
            All of the constants present in the equations of motion. The
            order does not matter.
        mass_matrix : sympy.Matrix, shape(n, n), optional
            This can be either the "full" mass matrix as in [2] or the
            "minimal" mass matrix as in [3]. The rows and columns must be
            ordered to match the order of the coordinates and speeds. In the
            case of the full mass matrix, the speeds should always be
            ordered before the speeds, i.e. x = [q, u].
        coordinate_derivatives : sympy.Matrix, shape(m, 1), optional
            If the "minimal" mass matrix, form [3], is supplied, then this
            column vector represents the right hand side of the kinematical
            differential equations.
        specifieds : sequence of SymPy Functions
            The specified exogenous inputs to the system. These should be
            functions of time and the order does not matter.
        linear_sys_solver : string or function
            Specify either `numpy` or `scipy` to use the linear solvers
            provided in each package or supply a function that solves a
            linear system Ax=b with the call signature x = solve(A, b). For
            example, if you need to use custom kwargs for the SciPy solver,
            pass in a lambda function that wraps the solver and sets them.
        constants_arg_type : string
            The generated function accepts two different types of arguments
            for the numerical values of the constants: either a ndarray of
            the constants values in the correct order or a dictionary
            mapping the constants symbols to the numerical values. If None,
            this is determined inside of the generated function and can
            cause a significant slow down for performance critical code. If
            you know apriori what arg types you need to support choose
            either ``array`` or ``dictionary``. Note that ``array`` is
            faster than ``dictionary``.
        specifieds_arg_type : string
            The generated function accepts three different types of
            arguments for the numerical values of the specifieds: either a
            ndarray of the specifieds values in the correct order, a
            function that generates the correctly ordered ndarray, or a
            dictionary mapping the specifieds symbols or tuples of thereof
            to floats, ndarrays, or functions. If None, this is determined
            inside of the generated function and can cause a significant
            slow down for performance critical code. If you know apriori
            what arg types you want to support choose either ``array``,
            ``function``, or ``dictionary``. The speed of each, from fast to
            slow, are ``array``, ``function``, ``dictionary``, None.

        """

        self.right_hand_side = right_hand_side
        self.coordinates = coordinates
        self.speeds = speeds
        self.constants = constants
        self.mass_matrix = mass_matrix
        self.coordinate_derivatives = coordinate_derivatives
        self.specifieds = specifieds
        self.linear_sys_solver = linear_sys_solver
        self.constants_arg_type = constants_arg_type
        self.specifieds_arg_type = specifieds_arg_type

        # As the order of the constants and specifieds arguments if not
        # important, allow Sets to be used as input. However, the order must be
        # maintained and converted to a Sequence.
        if constants is not None and not isinstance(constants, Sequence):
            self.constants = tuple(constants)
        if specifieds is not None and not isinstance(specifieds, Sequence):
            self.specifieds = tuple(specifieds)

        self.system_type = self._deduce_system_type(
            mass_matrix=mass_matrix,
            coordinate_derivatives=coordinate_derivatives)

        self.num_coordinates = len(coordinates)
        self.num_speeds = len(speeds)
        self.num_states = self.num_coordinates + self.num_speeds
        self.num_constants = len(constants)

        if self.specifieds is None:
            self.num_specifieds = 0
            self.specifieds_arg_type = None
        else:
            self.num_specifieds = len(specifieds)

        # These are pre-allocated storage for the numerical values used in
        # some of the rhs() evaluations.
        self._constants_values = np.empty(self.num_constants)
        self._specifieds_values = np.empty(self.num_specifieds)

        self._check_system_consitency()

    @property
    def linear_sys_solver(self):
        return self._linear_sys_solver

    @linear_sys_solver.setter
    def linear_sys_solver(self, v):

        if isinstance(v, type(lambda x: x)):
            self._solve_linear_system = v
            self._linear_sys_solver = v
        elif v == 'numpy':
            self._solve_linear_system = numpy.linalg.solve
            self._linear_sys_solver = v
        elif v == 'scipy':
            self._solve_linear_system = scipy.linalg.solve
            self._linear_sys_solver = v
        elif v == 'sympy':
            # dummy function
            self._solve_linear_system = lambda A, b: np.nan*np.ones_like(b)
            self._linear_sys_solver = v
        else:
            msg = '{} is not a valid solver.'
            raise ValueError(msg.format(self.linear_sys_solver))

    def _check_system_consitency(self):

        if self.system_type == 'min mass matrix':

            nr, nc = self.mass_matrix.shape
            assert self.num_speeds == nr == nc
            assert self.num_speeds == self.right_hand_side.shape[0]
            assert self.num_coordinates == self.coordinate_derivatives.shape[0]

        elif self.system_type == 'full mass matrix':

            nr, nc = self.mass_matrix.shape
            assert self.num_states == nr == nc
            assert self.num_states == self.right_hand_side.shape[0]
            assert self.coordinate_derivatives is None

        elif self.system_type == 'full rhs':

            assert self.num_states == self.right_hand_side.shape[0]
            assert self.mass_matrix is None
            assert self.coordinate_derivatives is None

    @staticmethod
    def list_syms(indent, syms):
        """Returns a string representation of a valid rst list of the
        symbols in the sequence syms and indents the list given the integer
        number of indentations."""
        indentation = ' ' * indent
        lst = '- ' + ('\n' + indentation + '- ').join([str(s) for s in syms])
        return indentation + lst

    def _convert_constants_dict_to_array(self, p):
        """Returns an array of numerical values from the constants
        dictionary in the correct order."""

        # NOTE : It's unfortunate that this has to be run at every rhs eval,
        # because subsequent calls to rhs() doesn't require different
        # constants. I suppose you can sub out all the constants in the EoMs
        # before passing them into the generator. That would beg for the
        # capability to support self.constants=None to skip all of this
        # stuff in the rhs eval.
        for i, c in enumerate(self.constants):
            self._constants_values[i] = p[c]

        return self._constants_values

    def _parse_constants(self, *args):
        """Returns an ndarray containing the numerical values of the
        constants in the correct order. If the constants are already an
        array, that array is returned."""

        p = args[-1]
        try:
            # NOTE : This emits "VisibleDeprecationWarning: using a non-integer
            # number instead of an integer will result in an error in the
            # future" along with the IndexError in NumPy 1.11. Not sure why it
            # gives the warning, it already gives and error with trying to
            # index with a SymPy symbol.
            p = self._convert_constants_dict_to_array(p)
        except IndexError:
            # p is an array so just return the args
            return args
        else:
            return args[:-1] + (p,)

    def _convert_specifieds_dict_to_array(self, x, t, r):

        for k, v in r.items():
            # TODO : Not sure if this is the best check here.
            if (isinstance(type(k), UndefinedFunction) or
                isinstance(k, Derivative)):
                k = (k,)
            idx = [self.specifieds.index(symmy) for symmy in k]
            try:
                self._specifieds_values[idx] = v(x, t)
            except TypeError:  # not callable
                # If not callable, then it should be a float, ndarray,
                # or indexable.
                self._specifieds_values[idx] = v

        return self._specifieds_values

    def _parse_specifieds(self, x, t, r, p):

        if isinstance(r, dict):
            # NOTE : This function sets self._specifieds_values, so here we
            # return nothing.
            self._convert_specifieds_dict_to_array(x, t, r)
        else:
            # More efficient.
            try:
                self._specifieds_values = r(x, t)
            except TypeError:  # not callable.
                # If not callable, then it should be a float or ndarray.
                self._specifieds_values = r

        return x, t, self._specifieds_values, p

    def _parse_all_args(self, *args):
        """Returns args formatted for the post 0.3.0 generators using all of
        the parsers. This is the slowest method and is used by default if no
        information is provided by the user on which type of args will be
        passed in."""

        args = self._parse_constants(*args)

        if self.specifieds is not None:
            args = self._parse_specifieds(*args)

        return args

    def _generate_rhs_docstring(self):

        template_values = {'num_states': self.num_states,
                           'state_list': self.list_syms(8, self.coordinates
                                                        + self.speeds),
                           'specified_call_sig': '',
                           'constants_explanation':
                               self._constants_doc_templates[
                                   self.constants_arg_type].format(**{
                                       'num_constants': self.num_constants,
                                       'constant_list': self.list_syms(
                                           8, self.constants)}),
                           'specifieds_explanation': ''}

        if self.specifieds is not None:
            template_values['specified_call_sig'] = ' r,'
            specified_template_values = {
                'num_specified': self.num_specifieds,
                'specified_list': self.list_syms(8, self.specifieds)}
            template_values['specifieds_explanation'] = \
                self._specifieds_doc_templates[self.specifieds_arg_type].format(
                    **specified_template_values)

        return self._rhs_doc_template.format(**template_values)

    def _create_rhs_function(self):
        """Returns a function in the form expected by scipy.integrate.odeint
        that computes the derivatives of the states."""

        p_arg_type = self.constants_arg_type
        r_arg_type = self.specifieds_arg_type

        if p_arg_type is None and r_arg_type is None:
            def rhs(*args):
                # args: x, t, p
                # or
                # args: x, t, r, p

                args = self._parse_all_args(*args)

                q = args[0][:self.num_coordinates]
                u = args[0][self.num_coordinates:]

                if self.constants:
                    xdot = self._base_rhs(q, u, *args[2:])
                else:
                    xdot = self._base_rhs(q, u, *(args[2:3] + ([],)))
                return xdot

            rhs.__doc__ = self._generate_rhs_docstring()

            return rhs

        if p_arg_type == 'dictionary':
            p = lambda *li : self._convert_constants_dict_to_array(li[-1])

        else:
            p = lambda *li : self._parse_constants(*li)[-1]

        if r_arg_type is None:
            if self.specifieds is not None:
                r = lambda *li : self._parse_specifieds(*li)[-2]

        elif r_arg_type == 'array':
            r = lambda *li : (li)[-2]

        elif r_arg_type == 'dictionary':
            r = lambda *li : self._convert_specifieds_dict_to_array(*li[:3])

        elif r_arg_type == 'function':
            r = lambda *li: li[2](*li[:2])

        def rhs(*args):
            # args: x, t, p
            # or
            # args: x, t, r, p

            q = args[0][:self.num_coordinates]
            u = args[0][self.num_coordinates:]

            if self.specifieds is None:
                if self.constants:
                    return self._base_rhs(q, u, p(*args))
                else:
                    return self._base_rhs(q, u, [])
            else:
                if self.constants:
                    return self._base_rhs(q, u, r(*args), p(*args))
                else:
                    return self._base_rhs(q, u, r(*args), [])

        rhs.__doc__ = self._generate_rhs_docstring()

        return rhs

    def _create_base_rhs_function(self):
        """Sets the self._base_rhs function. This function accepts arguments
        in this form: (q, u, p) or (q, u, r, p)."""

        if self.system_type == 'full rhs':

            self._base_rhs = self.eval_arrays

        elif (self.system_type == 'full mass matrix' and
              self.linear_sys_solver=='sympy'):

            def base_rhs(*args):
                M, xdot = self.eval_arrays(*args)
                return xdot

            self._base_rhs = base_rhs

        elif self.system_type == 'full mass matrix':

            def base_rhs(*args):

                M, F = self.eval_arrays(*args)
                return self._solve_linear_system(M, F)

            self._base_rhs = base_rhs

        elif (self.system_type == 'min mass matrix' and
              self.linear_sys_solver=='sympy'):

            xdot = np.empty(self.num_states, dtype=float)

            def base_rhs(*args):
                M, udot, qdot = self.eval_arrays(*args)
                xdot[:self.num_coordinates] = qdot
                xdot[self.num_coordinates:] = udot
                return xdot

            self._base_rhs = base_rhs

        elif self.system_type == 'min mass matrix':

            xdot = np.empty(self.num_states, dtype=float)

            def base_rhs(*args):
                M, F, qdot = self.eval_arrays(*args)
                if self.num_speeds == 1:
                    udot = F / M
                else:
                    udot = self._solve_linear_system(M, F)
                xdot[:self.num_coordinates] = qdot
                xdot[self.num_coordinates:] = udot
                return xdot

            self._base_rhs = base_rhs

    def define_inputs(self):
        """Sets self.inputs to the list of sequences [q, u, p] or [q, u, r,
        p]."""

        self.inputs = [self.coordinates, self.speeds, self.constants]
        if self.specifieds is not None:
            self.inputs.insert(2, self.specifieds)

    def generate(self):
        """Returns a function that evaluates the right hand side of the
        first order ordinary differential equations in one of two forms:

            x' = f(x, t, p)

            or

            x' = f(x, t, r, p)

        See the docstring of the generated function for more details.

        """

        if self.system_type == 'full rhs':
            self.generate_full_rhs_function()
        elif self.system_type == 'full mass matrix':
            self.generate_full_mass_matrix_function()
        elif self.system_type == 'min mass matrix':
            self.generate_min_mass_matrix_function()

        self._create_base_rhs_function()

        return self._create_rhs_function()


class CythonODEFunctionGenerator(ODEFunctionGenerator):

    def __init__(self, *args, **kwargs):

        self._options = {'tmp_dir': None,
                         'prefix': 'pydy_codegen',
                         'cse': True,
                         'verbose': False,
                         }
        for k, v in self._options.items():
            self._options[k] = kwargs.pop(k, v)

        if Cython is None:
            raise ImportError('Cython must be installed to use this class.')
        else:
            super(CythonODEFunctionGenerator, self).__init__(*args, **kwargs)

    __init__.__doc__ = ODEFunctionGenerator.__init__.__doc__

    def _cythonize(self, outputs, inputs):
        g = CythonMatrixGenerator(inputs, outputs,
                                  prefix=self._options['prefix'],
                                  cse=self._options['cse'])
        return g.compile(tmp_dir=self._options['tmp_dir'],
                         verbose=self._options['verbose'])

    def _cythonize_symbolic_lusolve(self, outputs, inputs):
        if not self._options['cse']:
            msg = 'cse has to be True if using the sympy linear system solver'
            raise ValueError(msg)
        g = CythonMatrixGenerator(inputs, outputs,
                                  prefix=self._options['prefix'],
                                  cse=True)
        # patch in the special generator
        g.c_matrix_generator = _CLUsolveGenerator(inputs, outputs)
        return g.compile(tmp_dir=self._options['tmp_dir'],
                         verbose=self._options['verbose'])

    def _set_eval_array(self, f):

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: f(q, u, p, *self._empties)
        else:
            self.eval_arrays = lambda q, u, r, p: f(q, u, r, p,
                                                    *self._empties)

    def generate_full_rhs_function(self):

        self.define_inputs()
        outputs = [self.right_hand_side]

        self._empties = (np.empty(self.num_states, dtype=float),)

        self._set_eval_array(self._cythonize(outputs, self.inputs))

    def generate_full_mass_matrix_function(self):

        self.define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side]

        mass_matrix_result = np.empty(self.num_states ** 2, dtype=float)
        rhs_result = np.empty(self.num_states, dtype=float)

        self._empties = (mass_matrix_result, rhs_result)

        if self.linear_sys_solver == 'sympy':
            self._set_eval_array(self._cythonize_symbolic_lusolve(outputs,
                                                                  self.inputs))
        else:
            self._set_eval_array(self._cythonize(outputs, self.inputs))

    def generate_min_mass_matrix_function(self):

        self.define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        mass_matrix_result = np.empty(self.num_speeds ** 2, dtype=float)
        rhs_result = np.empty(self.num_speeds, dtype=float)
        kin_diffs_result = np.empty(self.num_coordinates, dtype=float)
        self._empties = (mass_matrix_result, rhs_result, kin_diffs_result)

        if self.linear_sys_solver == 'sympy':
            self._set_eval_array(self._cythonize_symbolic_lusolve(outputs,
                                                                  self.inputs))
        else:
            self._set_eval_array(self._cythonize(outputs, self.inputs))


class LambdifyODEFunctionGenerator(ODEFunctionGenerator):

    def __init__(self, *args, **kwargs):

        self._options = {'cse': True}

        for k, v in self._options.items():
            self._options[k] = kwargs.pop(k, v)

        super(LambdifyODEFunctionGenerator, self).__init__(*args, **kwargs)

    __init__.__doc__ = ODEFunctionGenerator.__init__.__doc__

    def _lambdify(self, outputs):
        # TODO : We could forgo this substitution for generation speed
        # purposes and have lots of args for lambdify (like it used to be
        # done) but there may be some limitations on number of args.
        subs = {}
        vec_inputs = []
        if self.specifieds is None:
            def_vecs = ['q', 'u', 'p']
        else:
            def_vecs = ['q', 'u', 'r', 'p']

        for syms, vec_name in zip(self.inputs, def_vecs):
            v = sm.DeferredVector(vec_name)
            for i, sym in enumerate(syms):
                subs[sym] = v[i]
            vec_inputs.append(v)

        try:
            outputs = [me.msubs(output, subs) for output in outputs]
        except AttributeError:
            # msubs doesn't exist in SymPy < 0.7.6.
            outputs = [output.subs(subs) for output in outputs]

        modules = [{'ImmutableMatrix': np.array}, 'numpy']

        try:
            return sm.lambdify(vec_inputs, outputs, modules=modules,
                               cse=self._options['cse'])
        except TypeError:  # cse added in SymPy 1.9
            return sm.lambdify(vec_inputs, outputs, modules=modules)

    def generate_full_rhs_function(self):

        self.define_inputs()
        outputs = [self.right_hand_side]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: np.squeeze(f(q, u, p))
        else:
            self.eval_arrays = lambda q, u, r, p: np.squeeze(f(q, u, r, p))

    def generate_full_mass_matrix_function(self):

        self.define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: tuple([np.squeeze(o) for o in
                                                      f(q, u, p)])
        else:
            self.eval_arrays = lambda q, u, r, p: tuple([np.squeeze(o) for o
                                                         in f(q, u, r, p)])

    def generate_min_mass_matrix_function(self):

        self.define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: tuple([np.squeeze(o) for o in
                                                      f(q, u, p)])
        else:
            self.eval_arrays = lambda q, u, r, p: tuple([np.squeeze(o) for o
                                                         in f(q, u, r, p)])


class TheanoODEFunctionGenerator(ODEFunctionGenerator):

    def __init__(self, *args, **kwargs):

        if theano is None:
            raise ImportError('Theano must be installed to use this class.')
        else:
            super(TheanoODEFunctionGenerator, self).__init__(*args, **kwargs)

    __init__.__doc__ = ODEFunctionGenerator.__init__.__doc__

    def define_inputs(self):
        # Theano's input requires a flatted sequence instead of sequence of
        # sequences.
        specifieds = []
        if self.specifieds is not None:
            specifieds = self.specifieds
        self.inputs = chain(self.coordinates, self.speeds,
                            specifieds, self.constants)

    def _theanoize(self, outputs):

        self.define_inputs()

        old_check_input = theano.config.check_input
        old_allow_gc = theano.config.allow_gc
        try:
            # This affects compilation and removes the input check at each step.
            theano.config.check_input = False

            # Disable Theano garbage collection to lower the number of allocations.
            theano.config.allow_gc = False

            f_imp = theano_function(self.inputs, outputs,
                                    on_unused_input='ignore',
                                    mode=theano.Mode(linker='c'))
        finally:
            theano.config.check_input = old_check_input
            theano.config.allow_gc = old_allow_gc

        # While denoting an input as trusted lowers Theano overhead:
        #     f.trust_input = True
        # we can bypass additional overhead with the following function:
        def f(*args):
            for i in range(len(args)):
                f_imp.input_storage[i].storage[0] = args[i]
            f_imp.fn()
            return [f_imp.output_storage[i].data for i in range(len(outputs))]

        return f

    def generate_full_rhs_function(self):

        outputs = [self.right_hand_side]

        f = self._theanoize(outputs)

        def eval_arrays(*args):
            vals = map(np.asarray, np.hstack(args))
            return np.squeeze(f(*vals))

        self.eval_arrays = eval_arrays

    def generate_full_mass_matrix_function(self):

        outputs = [self.mass_matrix, self.right_hand_side]

        f = self._theanoize(outputs)

        def eval_arrays(*args):
            vals = map(np.asarray, np.hstack(args))
            return tuple([np.squeeze(o) for o in f(*vals)])

        self.eval_arrays = eval_arrays

    def generate_min_mass_matrix_function(self):

        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        f = self._theanoize(outputs)

        def eval_arrays(*args):
            vals = map(np.asarray, np.hstack(args))
            return tuple([np.squeeze(o) for o in f(*vals)])

        self.eval_arrays = eval_arrays


def generate_ode_function(*args, **kwargs):
    """This is a function wrapper to the above classes. The docstring is
    automatically generated below."""

    generators = {'lambdify': LambdifyODEFunctionGenerator,
                  'cython': CythonODEFunctionGenerator,
                  'theano': TheanoODEFunctionGenerator}

    generator = kwargs.pop('generator', 'lambdify')

    try:
        # See if user passed in a custom class.
        g = generator(*args, **kwargs)
    except TypeError:
        # See if user passed in a string.
        try:
            Generator = generators[generator]
            g = Generator(*args, **kwargs)
        except KeyError:
            msg = '{} is not a valid generator.'.format(generator)
            raise NotImplementedError(msg)
        else:
            return g.generate()
    else:
        return g.generate()


_docstr = ODEFunctionGenerator.__init__.__doc__
_extra_parameters_doc = \
"""\
        generator : string or and ODEFunctionGenerator, optional
            The method used for generating the numeric right hand side. The
            string options are {'lambdify'|'theano'|'cython'} with
            'lambdify' being the default. You can also pass in a custom
            subclass of ODEFunctionGenerator.

        Returns
        =======
        rhs : function
            A function which evaluates the derivaties of the states. See the
            function's docstring for more details after generation.
"""
generate_ode_function.__doc__ = ('' * 4 + _docstr + _extra_parameters_doc)
