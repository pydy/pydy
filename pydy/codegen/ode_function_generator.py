#!/usr/bin/env python

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

            q = args[0][:self.num_coordinates]
            u = args[0][self.num_coordinates:]

            if self.specifieds is not None:
                args = self._parse_specifieds(*args)

            xdot = self.eval_arrays(q, u, *args[2:])

            return xdot

        return rhs


class FullMassMatrixMixin(MixinBase):

    def create_rhs_function(self):

        def rhs(*args):
            # args: x, t, p
            # args: x, t, r, p

            q = args[0][:self.num_coordinates]
            u = args[0][self.num_coordinates:]

            if self.specifieds is not None:
                args = self._parse_specifieds(*args)

            M, F = self.eval_arrays(q, u, *args[2:])

            xdot = self.solve_linear_system(M, F)

            return xdot

        return rhs


class MinMassMatrixMixin(MixinBase):

    def create_rhs_function(self):

        xdot = np.empty(self.num_states, dtype=float)

        def rhs(*args):
            # args: x, t, p
            # args: x, t, r, p

            q = args[0][:self.num_coordinates]
            u = args[0][self.num_coordinates:]

            if self.specifieds is not None:
                args = self._parse_specifieds(*args)

            M, F, qdot = self.eval_arrays(q, u, *args[2:])

            if self.num_speeds == 1:
                udot = F / M
            else:
                udot = self.solve_linear_system(M, F)

            xdot[:self.num_coordinates] = qdot
            xdot[self.num_coordinates:] = udot

            return xdot

        return rhs


class ODEFunctionGenerator(object):

    @staticmethod
    def _deduce_system_type(**kwargs):

        if kwargs.pop('coordinate_derivatives') is not None:
            system_type = 'min mass matrix'
        elif kwargs.pop('mass_matrix') is not None:
            system_type = 'full mass matrix'
        else:
            system_type = 'full rhs'

        return system_type

    def __new__(cls,
                right_hand_side,
                coordinates,
                speeds,
                constants,
                mass_matrix=None,
                coordinate_derivatives=None,
                specifieds=None,
                linear_sys_solver='numpy'):
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

    def __init__(self,
                 right_hand_side,
                 coordinates,
                 speeds,
                 constants,
                 mass_matrix=None,
                 coordinate_derivatives=None,
                 specifieds=None,
                 linear_sys_solver='numpy'):
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
            linear system Ax=b with the call signature x = solve(A, b).

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
        # TODO : This fails for multiple calls if the tmp_dir is not set.
        # The module counter must not be advancing.
        return CythonMatrixGenerator(inputs, outputs).compile(tmp_dir='booger')

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
