#!/usr/bin/env python

import numpy as np
import scipy as sp
import sympy as sm

from ...models import multi_mass_spring_damper
from ..ode_function_generator import (ODEFunctionGenerator,
                                      CythonODEFunctionGenerator,
                                      LambdifyODEFunctionGenerator)


class TestODEFunctionGenerator(object):

    def setup(self):

        self.sys = multi_mass_spring_damper(2)

    def test_init(self):

        g = ODEFunctionGenerator(self.sys.eom_method.rhs(),
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols)

        assert g.num_coordinates == 2
        assert g.num_speeds == 2
        assert g.num_states == 4

    def test_init_full_rhs(self):

        g = ODEFunctionGenerator(self.sys.eom_method.rhs(),
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols)

        assert g.system_type == 'full rhs'

    def test_init_full_mass_matrix(self):

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full)

        assert g.system_type == 'full mass matrix'

    def test_init_min_mass_matrix(self):

        kin_diff_eqs = self.sys.eom_method.kindiffdict()
        coord_derivs = sm.Matrix([kin_diff_eqs[c.diff()] for c in
                                  self.sys.coordinates])

        g = ODEFunctionGenerator(self.sys.eom_method.forcing,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix,
                                 coordinate_derivatives=coord_derivs)

        assert g.system_type == 'min mass matrix'

    def test_set_linear_system_solver(self):

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full)

        assert g.solve_linear_system == np.linalg.solve

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver='numpy')

        assert g.solve_linear_system == np.linalg.solve

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver='scipy')

        assert g.solve_linear_system == sp.linalg.solve

        solver = lambda A, b: np.dot(np.inv(A),b)

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver=solver)

        assert g.solve_linear_system == solver


class TestODEFunctionGeneratorSubclasses(object):

    ode_function_subclasses = [CythonODEFunctionGenerator,
                               LambdifyODEFunctionGenerator]

    def setup(self):

        self.sys = multi_mass_spring_damper()

    def eval_rhs(self, rhs):

        xdot = rhs(np.array([1.0, 2.0]), 0.0, np.array([1.0, 2.0, 3.0]))

        expected_xdot = np.array([2.0, (-2.0 * 2.0 - 3.0 * 1.0) / 1.0])

        np.testing.assert_allclose(xdot, expected_xdot)

    def test_generate_full_rhs(self):

        rhs = self.sys.eom_method.rhs()

        for Subclass in self.ode_function_subclasses:

            g = Subclass(rhs,
                         self.sys.coordinates,
                         self.sys.speeds,
                         self.sys.constants_symbols)

            rhs_func = g.generate()

            self.eval_rhs(rhs_func)

    def test_generate_full_mass_matrix(self):

        for Subclass in self.ode_function_subclasses:

            g = Subclass(self.sys.eom_method.forcing_full,
                         self.sys.coordinates,
                         self.sys.speeds,
                         self.sys.constants_symbols,
                         mass_matrix=self.sys.eom_method.mass_matrix_full)

            rhs_func = g.generate()

            self.eval_rhs(rhs_func)

    def test_generate_min_mass_matrix(self):

        kin_diff_eqs = self.sys.eom_method.kindiffdict()
        coord_derivs = sm.Matrix([kin_diff_eqs[c.diff()] for c in
                                  self.sys.coordinates])

        for Subclass in self.ode_function_subclasses:

            g = Subclass(self.sys.eom_method.forcing,
                         self.sys.coordinates,
                         self.sys.speeds,
                         self.sys.constants_symbols,
                         mass_matrix=self.sys.eom_method.mass_matrix,
                         coordinate_derivatives=coord_derivs)

            rhs_func = g.generate()

            self.eval_rhs(rhs_func)
