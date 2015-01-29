#!/usr/bin/env python

import numpy as np
import numpy.linalg
import scipy.linalg
import sympy as sm
import sympy.physics.mechanics as me

from .cython_code import CythonMatrixGenerator

# Required to replace current code.
# TODO : Support specifieds!
# TODO : Create Generator class for lambdify.
# TODO : Create Generator class for Theano.

# Enhancements
# TODO : Support creating a function for the Jacobian of the rhs.
# TODO : Support functions of time in the ODEs.
# TODO : Create Generator class for Octave.
# TODO : Support output equations.


class FullRHSMixin(object):

    def create_rhs_function(self):
        """Returns a function in the form expected by scipy.integrate.odeint
        that computes the derivatives of the states."""

        def rhs(x, t, p):

            q = x[:self.num_coordinates]
            u = x[self.num_coordinates:]

            xdot = self.eval_arrays(q, u, p)

            return xdot

        return rhs


class FullMassMatrixMixin(object):

    def create_rhs_function(self):

        def rhs(x, t, p):

            q = x[:self.num_coordinates]
            u = x[self.num_coordinates:]

            M, F = self.eval_arrays(q, u, p)

            xdot = self.solve_linear_system(M, F)

            return xdot

        return rhs


class MinMassMatrixMixin(object):

    def create_rhs_function(self):

        xdot = np.empty(self.num_states, dtype=float)

        def rhs(x, t, p):

            q = x[:self.num_coordinates]
            u = x[self.num_coordinates:]
            M, F, qdot = self.eval_arrays(q, u, p)
            if F.shape is tuple() or F.shape[0]:
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
        specifieds : ssequence of SymPy Functions
            The specified exogneous inputs to the ODEs.
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

        def f():
            pass

        if isinstance(self.linear_sys_solver, type(f)):
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

    def generate(self):
        """Returns a function that evaluates the right hand side of the
        first order odes.

        x' = f(x, t, p)

        """
        if self.system_type == 'full rhs':
            self.generate_full_rhs_function()
        elif self.system_type == 'full mass matrix':
            self.generate_full_mass_matrix_function()
        elif self.system_type == 'min mass matrix':
            self.generate_min_mass_matrix_function()

        return self.create_rhs_function()


class CythonODEFunctionGenerator(ODEFunctionGenerator):

    @staticmethod
    def _cythonize(outputs, inputs):
        # TODO : This fails for multiple calls if the tmp_dir is not set.
        # The module counter must not be advancing.
        return CythonMatrixGenerator(outputs, inputs).compile(tmp_dir='booger')

    def generate_full_rhs_function(self):

        outputs = [self.right_hand_side]
        inputs = [self.coordinates, self.speeds, self.constants]

        rhs_result = np.empty(self.num_states, dtype=float)

        f = self._cythonize(outputs, inputs)

        self.eval_arrays = lambda q, u, p: f(q, u, p, rhs_result)

    def generate_full_mass_matrix_function(self):

        outputs = [self.mass_matrix, self.right_hand_side]
        inputs = [self.coordinates, self.speeds, self.constants]

        f = self._cythonize(outputs, inputs)

        mass_matrix_result = np.empty(self.num_states ** 2, dtype=float)
        rhs_result = np.empty(self.num_states, dtype=float)

        self.eval_arrays = lambda q, u, p: f(q, u, p, mass_matrix_result,
                                             rhs_result)

    def generate_min_mass_matrix_function(self):

        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]
        inputs = [self.coordinates, self.speeds, self.constants]

        f = self._cythonize(outputs, inputs)

        mass_matrix_result = np.empty(self.num_speeds ** 2, dtype=float)
        rhs_result = np.empty(self.num_speeds, dtype=float)
        kin_diffs_result = np.empty(self.num_coordinates, dtype=float)

        self.eval_arrays = lambda q, u, p: f(q, u, p, mass_matrix_result,
                                             rhs_result, kin_diffs_result)


class LambdifyODEFunctionGenerator(ODEFunctionGenerator):

    @staticmethod
    def _lambdify(inputs, outputs):
        # TODO : We could forgo this substitution for speed purposes and
        # have lots of args for lambdify (like it used to be done) but there
        # may be some limitiations on number of args.
        subs = {}
        vec_inputs = []
        for syms, vec_name in zip(inputs, ['q', 'u', 'p']):
            v = sm.DeferredVector(vec_name)
            for i, sym in enumerate(syms):
                subs[sym] = v[i]
            vec_inputs.append(v)
        outputs = [me.msubs(output, subs) for output in outputs]

        modules = [{'ImmutableMatrix': np.array}, 'numpy']

        return sm.lambdify(vec_inputs, outputs, modules=modules)

    def generate_full_rhs_function(self):

        inputs = [self.coordinates, self.speeds, self.constants]
        outputs = [self.right_hand_side]

        f = self._lambdify(inputs, outputs)
        self.eval_arrays = lambda q, u, p: np.squeeze(f(q, u, p))

    def generate_full_mass_matrix_function(self):

        inputs = [self.coordinates, self.speeds, self.constants]
        outputs = [self.mass_matrix, self.right_hand_side]

        f = self._lambdify(inputs, outputs)
        self.eval_arrays = lambda q, u, p: tuple([np.squeeze(o) for o in
                                                  f(q, u, p)])

    def generate_min_mass_matrix_function(self):

        inputs = [self.coordinates, self.speeds, self.constants]
        outputs = [self.mass_matrix, self.right_hand_side,
                   self.coordinate_derivatives]

        f = self._lambdify(inputs, outputs)
        self.eval_arrays = lambda q, u, p: tuple([np.squeeze(o) for o in
                                                  f(q, u, p)])
