#!/usr/bin/env python

import numpy as np
import scipy as sp
import sympy as sm

Cython = sm.external.import_module('Cython')
theano = sm.external.import_module('theano')

from ... import models
from ..ode_function_generator import (ODEFunctionGenerator,
                                      LambdifyODEFunctionGenerator,
                                      CythonODEFunctionGenerator,
                                      TheanoODEFunctionGenerator)


class TestODEFunctionGenerator(object):

    def setup(self):

        self.sys = models.multi_mass_spring_damper(2)

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

        solver = lambda A, b: np.dot(np.inv(A), b)

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver=solver)

        assert g.solve_linear_system == solver


class TestODEFunctionGeneratorSubclasses(object):

    ode_function_subclasses = [LambdifyODEFunctionGenerator,
                               CythonODEFunctionGenerator,
                               TheanoODEFunctionGenerator]

    def setup(self):

        self.sys = models.multi_mass_spring_damper()

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

    def test_generate_with_specifieds(self):

        # This model has three specifieds.
        sys = models.n_link_pendulum_on_cart(2, True, True)

        kin_diff_eqs = sys.eom_method.kindiffdict()
        coord_derivs = sm.Matrix([kin_diff_eqs[c.diff()] for c in
                                  sys.coordinates])

        rhs = sys.eom_method.rhs()

        common_args = [sys.coordinates, sys.speeds, sys.constants_symbols]
        args_dct = {}
        args_dct['full rhs'] = [rhs] + common_args
        args_dct['full mass matrix'] = [sys.eom_method.forcing_full] + common_args
        args_dct['min mass matrix'] = [sys.eom_method.forcing] + common_args

        kwargs_dct = {}
        kwargs_dct['full rhs'] = {}
        kwargs_dct['full mass matrix'] = \
            {'mass_matrix': sys.eom_method.mass_matrix_full}
        kwargs_dct['min mass matrix'] = \
            {'mass_matrix': sys.eom_method.mass_matrix,
             'coordinate_derivatives': coord_derivs}

        for Subclass in self.ode_function_subclasses:
            for system_type, args in args_dct.items():

                g = Subclass(*args, specifieds=sys.specifieds_symbols,
                             **kwargs_dct[system_type])

                f = g.generate()

                rand = np.random.random

                x = rand(g.num_states)
                r = rand(g.num_specifieds)
                p = rand(g.num_constants)

                subs = {}
                for arr, syms in zip([x, r, p], [sys.states,
                                                 sys.specifieds_symbols,
                                                 sys.constants_symbols]):
                    for val, sym in zip(arr, syms):
                        subs[sym] = val

                expected = sm.matrix2numpy(rhs.subs(subs), dtype=float).squeeze()

                xdot = f(x, 0.0, r, p)

                np.testing.assert_allclose(xdot, expected)

    def test_rhs_args(self):

        # There are four specified inputs available.
        sys = models.n_link_pendulum_on_cart(3, True, True)
        right_hand_side = sys.eom_method.rhs()

        for Subclass in self.ode_function_subclasses:

            g = Subclass(right_hand_side,
                         sys.coordinates, sys.speeds,
                         sys.constants_symbols,
                         specifieds=sys.specifieds_symbols)

            rhs = g.generate()

            x = np.random.random(g.num_states)
            p = np.random.random(g.num_constants)

            r = dict(zip(g.specifieds, [1.0, 2.0, 3.0, 4.0]))

            xd_01 = rhs(x, 0.0, r, p)

            r = {tuple(g.specifieds):
                 lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])}
            xd_02 = rhs(x, 0.0, r, p)
            np.testing.assert_allclose(xd_01, xd_02)

            r = {g.specifieds[0]: lambda x, t: np.ones(1),
                 (g.specifieds[3], g.specifieds[1]):
                 lambda x, t: np.array([4.0, 2.0]),
                 g.specifieds[2]: 3.0 * np.ones(1)}
            xd_03 = rhs(x, 0.0, r, p)
            np.testing.assert_allclose(xd_01, xd_03)

            r = np.array([1.0, 2.0, 3.0, 4.0])
            xd_04 = rhs(x, 0.0, r, p)
            np.testing.assert_allclose(xd_01, xd_04)

            r = lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])
            xd_05 = rhs(x, 0.0, r, p)
            np.testing.assert_allclose(xd_01, xd_05)
