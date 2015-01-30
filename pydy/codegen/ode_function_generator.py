#!/usr/bin/env python

import warnings

import numpy as np
import numpy.linalg
import scipy.linalg
import sympy as sm
import sympy.physics.mechanics as me
from sympy.core.function import UndefinedFunction

Cython = sm.external.import_module('Cython')
theano = sm.external.import_module('theano')
if theano:
    from sympy.printing.theanocode import theano_function

from .cython_code import CythonMatrixGenerator


class MixinBase(object):

    def _parse_constants(self, *args):

        p = args[-1]
        try:
            p = np.array([p[c] for c in self.constants])
            return args[:-1] + (p,)
        except IndexError:
            return args

    def _parse_old_style_extra_args(self, *args):

        # DEPRECATED : Remove before 0.4.0 release.
        try:
            # old style always has three args: x, t, args and the last one
            # is a dictionary which at least contains the key 'constants'.
            # So if the last arg is so, then extract.
            args[-1]['constants']
        except (KeyError, IndexError):
            return args
        else:
            new_args = list(args[:-1])  # gets x and t

            if self.specifieds is not None:
                new_args.append(args[-1]['specified'])

            new_args.append(args[-1]['constants'])

            return tuple(new_args)

    def _parse_specifieds(self, x, t, r, p):

        if isinstance(r, dict):
            # TODO : This should be instatiated outside of the RHS function.
            specified_values = np.zeros(self.num_specifieds)

            for k, v in r.items():
                # TODO : Not sure if this is the best check here.
                if isinstance(type(k), UndefinedFunction):
                    k = (k,)
                idx = [self.specifieds.index(symmy) for symmy in k]
                try:
                    specified_values[idx] = v(x, t)
                except TypeError:  # not callable
                    # If not callable, then it should be a float, ndarray,
                    # or indexable.
                    specified_values[idx] = v
        else:
            # More efficient.
            try:
                specified_values = r(x, t)
            except TypeError:  # not callable.
                # If not callable, then it should be a float or ndarray.
                specified_values = r
            # If the value is just a float, then convert to a 1D array.
            try:
                len(specified_values)
            except TypeError:
                specified_values = np.asarray([specified_values])

        return x, t, specified_values, p


class FullRHSMixin(MixinBase):

    def create_rhs_function(self):
        """Returns a function in the form expected by scipy.integrate.odeint
        that computes the derivatives of the states."""

        def rhs(*args):
            # args: x, t, p
            # args: x, t, r, p

            args = self._parse_old_style_extra_args(*args)

            args = self._parse_constants(*args)

            if self.specifieds is not None:
                args = self._parse_specifieds(*args)

            q = args[0][:self.num_coordinates]
            u = args[0][self.num_coordinates:]

            xdot = self.eval_arrays(q, u, *args[2:])

            return xdot

        rhs.__doc__ = self._generate_rhs_docstring()

        return rhs


class FullMassMatrixMixin(MixinBase):

    def create_rhs_function(self):

        def rhs(*args):
            # args: x, t, p
            # args: x, t, r, p

            args = self._parse_old_style_extra_args(*args)

            args = self._parse_constants(*args)

            if self.specifieds is not None:
                args = self._parse_specifieds(*args)

            q = args[0][:self.num_coordinates]
            u = args[0][self.num_coordinates:]

            M, F = self.eval_arrays(q, u, *args[2:])

            xdot = self.solve_linear_system(M, F)

            return xdot

        rhs.__doc__ = self._generate_rhs_docstring()

        return rhs


class MinMassMatrixMixin(MixinBase):

    def create_rhs_function(self):

        xdot = np.empty(self.num_states, dtype=float)

        def rhs(*args):
            # args: x, t, p
            # args: x, t, r, p

            args = self._parse_old_style_extra_args(*args)

            args = self._parse_constants(*args)

            if self.specifieds is not None:
                args = self._parse_specifieds(*args)

            q = args[0][:self.num_coordinates]
            u = args[0][self.num_coordinates:]

            M, F, qdot = self.eval_arrays(q, u, *args[2:])

            if self.num_speeds == 1:
                udot = F / M
            else:
                udot = self.solve_linear_system(M, F)

            xdot[:self.num_coordinates] = qdot
            xdot[self.num_coordinates:] = udot

            return xdot

        rhs.__doc__ = self._generate_rhs_docstring()

        return rhs


class ODEFunctionGenerator(object):

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
    The current time.{specified_explanation}
p : dictionary len({num_constants}) or ndarray shape({num_constants},)
    The dictionary maps the constant symbols to floats. The constants are,
    in order:
{constant_list}

Returns
=======
dx : ndarray, shape({num_states},)
    The derivative of the state vector.

"""

    _specifieds_doc_template = \
"""
r : dictionary; ndarray, shape({num_specified},); function

    There are three options for this dictionary. (1) is more flexible but
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

    @staticmethod
    def _deduce_system_type(**kwargs):

        if kwargs.pop('coordinate_derivatives') is not None:
            system_type = 'min mass matrix'
        elif kwargs.pop('mass_matrix') is not None:
            system_type = 'full mass matrix'
        else:
            system_type = 'full rhs'

        return system_type

    def __new__(cls, right_hand_side, coordinates, speeds, constants,
                mass_matrix=None, coordinate_derivatives=None,
                specifieds=None, linear_sys_solver='numpy'):
        """Returns an instance of the class with the appropriate mixin class
        based on what type of system was provided."""

        system_type = cls._deduce_system_type(
            mass_matrix=mass_matrix,
            coordinate_derivatives=coordinate_derivatives)

        if system_type == 'min mass matrix':
            bases = (MinMassMatrixMixin, cls)
        elif system_type == 'full mass matrix':
            bases = (FullMassMatrixMixin, cls)
        elif system_type == 'full rhs':
            bases = (FullRHSMixin, cls)

        return object.__new__(type(cls.__name__, bases, dict(cls.__dict__)))

    def __init__(self, right_hand_side, coordinates, speeds, constants,
                 mass_matrix=None, coordinate_derivatives=None,
                 specifieds=None, linear_sys_solver='numpy'):
        """

        Parameters
        ==========
        right_hand_side : SymPy Matrix, shape(n, 1)
            A column vector containing the symbolic expressions for the
            right hand side of the ordinary differential equations. If the
            right hand side has been solved for symbolically then x' = f(x,
            t, r, p), if it hasn't then you must supply the accompanying
            mass matrix M(x) * x' = f(x, t, r, p).
        coordinates : sequence of SymPy Functions
            The generalized coordinates.
        speeds : sequence of SymPy Functions
            The generalized speeds.
        constants : sequence of SymPy Symbols
            All of the constants present in the right hand side and mass
            matrix.
        mass_matrix : sympy.Matrix, shape(n, n), optional
            This can be either the "full" mass matrix in M(x) * x' = f(x, t,
            r, p) or the minimal mass matrix in M(q) * u' = f(u, q, t, r,
            p).
        coordinate_derivatives : sympy.Matrix, shape(m, 1), optional
            If the "minimal" mass matrix is supplied, then this matrix
            represents the right hand side of q' = g(q, t, p).
        specifieds : sequence of SymPy Functions
            The specified exogenous inputs to the ODEs. These should be
            functions of time.
        linear_sys_solver : string or function
            Specify either `numpy` or `scipy` to use the linear solvers
            provided in each package or supply a function that solves a
            linear system Ax=b with the call signature x = solve(A, b). If
            you need to use custom kwargs for the scipy solver, pass in a
            lambda function that wraps the solver and sets them.

        """

        self.right_hand_side = right_hand_side
        self.coordinates = coordinates
        self.speeds = speeds
        self.constants = constants
        self.mass_matrix = mass_matrix
        self.coordinate_derivatives = coordinate_derivatives
        self.specifieds = specifieds
        self.linear_sys_solver = linear_sys_solver

        self.system_type = self._deduce_system_type(
            mass_matrix=mass_matrix,
            coordinate_derivatives=coordinate_derivatives)

        self.num_coordinates = len(coordinates)
        self.num_speeds = len(speeds)
        self.num_states = self.num_coordinates + self.num_speeds
        self.num_constants = len(constants)

        if self.specifieds is not None:
            self.num_specifieds = len(specifieds)
        else:
            self.num_specifieds = 0

        self._check_system_consitency()
        self._set_linear_sys_solver()

    def _set_linear_sys_solver(self):

        if isinstance(self.linear_sys_solver, type(lambda x: x)):
            self.solve_linear_system = self.linear_sys_solver
        elif self.linear_sys_solver == 'numpy':
            self.solve_linear_system = numpy.linalg.solve
        elif self.linear_sys_solver == 'scipy':
            self.solve_linear_system = scipy.linalg.solve
        else:
            msg = '{} is not a valid solver'
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
        indentation = ' ' * indent
        lst = '- ' + ('\n' + indentation + '- ').join([str(s) for s in syms])
        return indentation + lst

    def _generate_rhs_docstring(self):

        template_values = {'num_states': self.num_states,
                           'state_list': self.list_syms(8, self.coordinates
                                                        + self.speeds),
                           'num_constants': self.num_constants,
                           'constant_list': self.list_syms(8, self.constants),
                           'specified_call_sig': '',
                           'specified_explanation': ''}

        if self.specifieds is not None:
            template_values['specified_call_sig'] = ' r,'
            specified_template_values = {
                'num_specified': self.num_specifieds,
                'specified_list': self.list_syms(8, self.specifieds)}
            template_values['specified_explanation'] = \
                self._specifieds_doc_template.format(**specified_template_values)

        return self._rhs_doc_template.format(**template_values)

    def _define_inputs(self):

        self.inputs = [self.coordinates, self.speeds, self.constants]

        if self.specifieds is not None:
            self.inputs.insert(2, self.specifieds)

    def generate(self):
        """Returns a function that evaluates the right hand side of the
        first order ordinary differential equations in one of two forms:

        x' = f(x, t, p)

        or

        x' = f(x, t, r, p)

        """
        if self.system_type == 'full rhs':
            self.generate_full_rhs_function()
        elif self.system_type == 'full mass matrix':
            self.generate_full_mass_matrix_function()
        elif self.system_type == 'min mass matrix':
            self.generate_min_mass_matrix_function()

        return self.create_rhs_function()


class CythonODEFunctionGenerator(ODEFunctionGenerator):

    def __init__(self, *args, **kwargs):

        if Cython is None:
            raise Exception('Cython must be installed to use this class.')
        else:
            super(CythonODEFunctionGenerator, self).__init__(*args, **kwargs)

    @staticmethod
    def _cythonize(outputs, inputs):
        return CythonMatrixGenerator(inputs, outputs).compile()

    def _set_eval_array(self, f):

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: f(q, u, p, *self.empties)
        else:
            self.eval_arrays = lambda q, u, r, p: f(q, u, r, p, *self.empties)

    def generate_full_rhs_function(self):

        self._define_inputs()
        outputs = [self.right_hand_side]

        self.empties = (np.empty(self.num_states, dtype=float),)

        self._set_eval_array(self._cythonize(outputs, self.inputs))

    def generate_full_mass_matrix_function(self):

        self._define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side]

        mass_matrix_result = np.empty(self.num_states ** 2, dtype=float)
        rhs_result = np.empty(self.num_states, dtype=float)

        self.empties = (mass_matrix_result, rhs_result)

        self._set_eval_array(self._cythonize(outputs, self.inputs))

    def generate_min_mass_matrix_function(self):

        self._define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        mass_matrix_result = np.empty(self.num_speeds ** 2, dtype=float)
        rhs_result = np.empty(self.num_speeds, dtype=float)
        kin_diffs_result = np.empty(self.num_coordinates, dtype=float)
        self.empties = (mass_matrix_result, rhs_result, kin_diffs_result)

        self._set_eval_array(self._cythonize(outputs, self.inputs))


class LambdifyODEFunctionGenerator(ODEFunctionGenerator):

    def _lambdify(self, outputs):
        # TODO : We could forgo this substitution for speed purposes and
        # have lots of args for lambdify (like it used to be done) but there
        # may be some limitations on number of args.
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
        outputs = [me.msubs(output, subs) for output in outputs]

        modules = [{'ImmutableMatrix': np.array}, 'numpy']

        return sm.lambdify(vec_inputs, outputs, modules=modules)

    def generate_full_rhs_function(self):

        self._define_inputs()
        outputs = [self.right_hand_side]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: np.squeeze(f(q, u, p))
        else:
            self.eval_arrays = lambda q, u, r, p: np.squeeze(f(q, u, r, p))

    def generate_full_mass_matrix_function(self):

        self._define_inputs()
        outputs = [self.mass_matrix, self.right_hand_side]

        f = self._lambdify(outputs)

        if self.specifieds is None:
            self.eval_arrays = lambda q, u, p: tuple([np.squeeze(o) for o in
                                                      f(q, u, p)])
        else:
            self.eval_arrays = lambda q, u, r, p: tuple([np.squeeze(o) for o
                                                         in f(q, u, r, p)])

    def generate_min_mass_matrix_function(self):

        self._define_inputs()
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
            raise Exception('Theano must be installed to use this class.')
        else:
            super(TheanoODEFunctionGenerator, self).__init__(*args, **kwargs)

    def _define_inputs(self):

        if self.specifieds is None:
            self.inputs = self.coordinates + self.speeds + self.constants
        else:
            self.inputs = (self.coordinates + self.speeds + self.specifieds
                           + self.constants)

    def _theanoize(self, outputs):

        self._define_inputs()

        f = theano_function(self.inputs, outputs, on_unused_input='ignore')

        # Theano will run faster if you trust the input. I'm not sure
        # what the implications of this are. See:
        # http://deeplearning.net/software/theano/tutorial/faq.html#faster-small-theano-function
        # Note that map(np.asarray, np.hstack(args)) is required if
        # trust_input is True. If it is False, then it will sanitize the
        # inputs. I'm not sure which one is faster.
        f.trust_input = True

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
    """Returns a numerical function which can evaluate the right hand side
    of the first order ordinary differential equations from a system
    described by one of the following three forms:

        [1] x' = F(x, t, r, p)

        [2] M(x, p) x' = F(x, t, r, p)

        [3] M(q, p) u' = F(q, u, t, r, p)
            q' = G(q, u, t, r, p)

    where

        x : states
        t : time
        r : specified inputs
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
        right hand side has been solved for symbolically then x' = f(x,
        t, r, p), if it hasn't then you must supply the accompanying
        mass matrix M(x) * x' = f(x, t, r, p).
    coordinates : sequence of SymPy Functions
        The generalized coordinates.
    speeds : sequence of SymPy Functions
        The generalized speeds.
    constants : sequence of SymPy Symbols
        All of the constants present in the right hand side and mass
        matrix.
    mass_matrix : sympy.Matrix, shape(n, n), optional
        This can be either the "full" mass matrix in M(x) * x' = f(x, t,
        r, p) or the minimal mass matrix in M(q) * u' = f(u, q, t, r,
        p).
    coordinate_derivatives : sympy.Matrix, shape(m, 1), optional
        If the "minimal" mass matrix is supplied, then this matrix
        represents the right hand side of q' = g(q, t, p).
    specifieds : sequence of SymPy Functions
        The specified exogenous inputs to the ODEs. These should be
        functions of time.
    linear_sys_solver : string or function
        Specify either `numpy` or `scipy` to use the linear solvers
        provided in each package or supply a function that solves a
        linear system Ax=b with the call signature x = solve(A, b). If
        you need to use custom kwargs for the scipy solver, pass in a
        lambda function that wraps the solver and sets them.
    generator : string, {'lambdify'|'theano'|'cython'}, optional
        The method used for generating the numeric right hand side.

    Returns
    -------
    rhs_func : function
        A function which evaluates the derivaties of the states. See the
        function's docstring for more details after generation.

    Notes
    -----

    This function also supports the pre-0.3.0 arguments for backwards
    compatibility which was strictly form [2] above:

    mass_matrix : sympy.Matrix, shape(n,n)
        The symbolic mass matrix of the system. The rows should correspond
        to the coordinates and speeds.
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


    """
    if len(args) == 5:
        warnings.warn("The pre-0.3.0 arguments for generate_ode_function are "
                      "deprecated. Please check the documentation for the new "
                      "style arguments.", DeprecationWarning)
        # These are the old style args, so we rearrange into the new style
        # args for backwards compatibility.
        mass_matrix, forcing_vector, constants, coordinates, speeds = args
        args = (forcing_vector, coordinates, speeds, constants)
        kwargs['mass_matrix'] = mass_matrix
        kwargs['specifieds'] = kwargs.pop('specified')

    generators = {'lambdify': LambdifyODEFunctionGenerator,
                  'cython': CythonODEGenerator,
                  'theano': TheanoODEGenerator}

    generator_type = kwargs.pop('generator')
    try:
        g = generators[generator_type](*args, **kwargs)
    except KeyError:
        msg = '{} is not a valid generator.'.format(generator_type)
        raise ValueError(msg)

    return g.generate()
