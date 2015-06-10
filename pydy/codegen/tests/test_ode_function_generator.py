#!/usr/bin/env python

from random import choice

import numpy as np
import scipy as sp
import sympy as sm

Cython = sm.external.import_module('Cython')
theano = sm.external.import_module('theano')

from ... import models
from ..ode_function_generators import (ODEFunctionGenerator,
                                       LambdifyODEFunctionGenerator,
                                       CythonODEFunctionGenerator,
                                       TheanoODEFunctionGenerator)


class TestODEFunctionGenerator(object):

    def setup(self):

        self.sys = models.multi_mass_spring_damper(2)
        self.rhs = self.sys.eom_method.rhs()
        self.generator = ODEFunctionGenerator(self.rhs,
                                              self.sys.coordinates,
                                              self.sys.speeds,
                                              self.sys.constants_symbols)

    def test_init_full_rhs(self):

        g = ODEFunctionGenerator(self.rhs,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols)

        assert g.num_coordinates == 2
        assert g.num_speeds == 2
        assert g.num_states == 4
        assert g.system_type == 'full rhs'


    def test_init_full_mass_matrix(self):

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full)

        assert g.num_coordinates == 2
        assert g.num_speeds == 2
        assert g.num_states == 4
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

        assert g.num_coordinates == 2
        assert g.num_speeds == 2
        assert g.num_states == 4
        assert g.system_type == 'min mass matrix'

    def test_set_linear_system_solver(self):

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full)

        assert g._solve_linear_system == np.linalg.solve

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver='numpy')

        assert g._solve_linear_system == np.linalg.solve

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver='scipy')

        assert g._solve_linear_system == sp.linalg.solve

        solver = lambda A, b: np.dot(np.inv(A), b)

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver=solver)

        assert g._solve_linear_system == solver


class TestODEFunctionGeneratorSubclasses(object):

    ode_function_subclasses = [LambdifyODEFunctionGenerator,
                               CythonODEFunctionGenerator,
                               TheanoODEFunctionGenerator]

    def setup(self):

        self.sys = models.multi_mass_spring_damper()
        # Best keep these in order, otherwise it may change between SymPy
        # versions.
        self.constants = list(sm.ordered(self.sys.constants_symbols))

    def eval_rhs(self, rhs):

        # In order:
        c0 = 1.0
        k0 = 2.0
        m0 = 3.0

        x0 = 1.0
        v0 = 2.0

        xdot = rhs(np.array([x0, v0]), 0.0, np.array([c0, k0, m0]))

        expected_xdot = np.array([v0, (-c0 * v0 - k0 * x0) / m0])

        np.testing.assert_allclose(xdot, expected_xdot)

    def test_init_doc(self):

        for Subclass in self.ode_function_subclasses:
            assert (Subclass.__init__.__doc__ ==
                    ODEFunctionGenerator.__init__.__doc__)

    def test_generate_full_rhs(self):

        rhs = self.sys.eom_method.rhs()

        for Subclass in self.ode_function_subclasses:

            g = Subclass(rhs,
                         self.sys.coordinates,
                         self.sys.speeds,
                         self.constants)

            rhs_func = g.generate()

            self.eval_rhs(rhs_func)

    def test_generate_full_mass_matrix(self):

        for Subclass in self.ode_function_subclasses:

            g = Subclass(self.sys.eom_method.forcing_full,
                         self.sys.coordinates,
                         self.sys.speeds,
                         self.constants,
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
                         self.constants,
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

                try:
                    expected = sm.matrix2numpy(rhs.subs(subs),
                                               dtype=float).squeeze()
                except TypeError:
                    # Earlier SymPy versions don't support the dtype kwarg.
                    expected = np.asarray(sm.matrix2numpy(rhs.subs(subs)),
                                          dtype=float).squeeze()

                xdot = f(x, 0.0, r, p)

                np.testing.assert_allclose(xdot, expected)

    def test_old_rhs_args(self):

        # DEPRECATED : This should be removed before the 0.4.0 release.

        # There are four specified inputs available.
        sys = models.n_link_pendulum_on_cart(3, True, True)
        right_hand_side = sys.eom_method.rhs()

        # It is only necessary to check one Generator because all of the
        # arg handling is shared among them.
        g = LambdifyODEFunctionGenerator(right_hand_side, sys.coordinates,
                                         sys.speeds, sys.constants_symbols,
                                         specifieds=sys.specifieds_symbols)

        rhs = g.generate()

        x = np.random.random(g.num_states)
        t = 0.0
        r = np.random.random(g.num_specifieds)
        p = np.random.random(g.num_constants)

        # Compute with new style args.
        xd_01 = rhs(x, t, r, p)

        # Now check to see if old style extra args works.
        args = {'specified': r, 'constants': p}

        xd_02 = rhs(x, t, args)
        np.testing.assert_allclose(xd_01, xd_02)

    def test_rhs_args(self):
        # This test takes a while to run but it checks all the combinations.

        # There are eight constants and four specified inputs available.
        sys = models.n_link_pendulum_on_cart(3, True, True)
        right_hand_side = sys.eom_method.rhs()
        constants = list(sm.ordered(sys.constants_symbols))
        specifieds = list(sm.ordered(sys.specifieds_symbols))

        constants_arg_types = [None, 'array', 'dictionary']
        specifieds_arg_types = [None, 'array', 'function', 'dictionary']

        p_array = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
        p_dct = dict(zip(constants, p_array))

        p = {}
        p[None] = choice([p_array, p_dct])
        p['array'] = p_array
        p['dictionary'] = p_dct

        r_array = np.array([1.0, 2.0, 3.0, 4.0])
        r_dct_1 = dict(zip(specifieds, r_array))
        r_dct_2 = {tuple(specifieds):
                   lambda x, t: r_array}
        r_dct_3 = {specifieds[0]: lambda x, t: np.ones(1),
                   (specifieds[3], specifieds[1]):
                   lambda x, t: np.array([4.0, 2.0]),
                   specifieds[2]: 3.0 * np.ones(1)}
        r_func = lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])

        r = {}
        r[None] = choice([r_array, r_dct_1, r_dct_2, r_dct_3, r_func])
        r['array'] = r_array
        r['dictionary'] = choice([r_dct_1, r_dct_2, r_dct_3])
        r['function'] = r_func

        x = np.random.random(len(sys.states))

        for p_arg_type in constants_arg_types:
            for r_arg_type in specifieds_arg_types:

                g = LambdifyODEFunctionGenerator(right_hand_side,
                                                 sys.coordinates,
                                                 sys.speeds,
                                                 constants,
                                                 specifieds=specifieds,
                                                 constants_arg_type=p_arg_type,
                                                 specifieds_arg_type=r_arg_type)
                rhs = g.generate()

                xdot = rhs(x, 0.0, r[r_arg_type], p[p_arg_type])

                try:
                    np.testing.assert_allclose(xdot, last_xdot)
                except NameError:
                    pass

                last_xdot = xdot

        # Now make sure it all works with specifieds=None
        sys = models.n_link_pendulum_on_cart(3, False, False)
        right_hand_side = sys.eom_method.rhs()
        constants = list(sm.ordered(sys.constants_symbols))

        del last_xdot

        for p_arg_type in constants_arg_types:
            for r_arg_type in specifieds_arg_types:

                g = LambdifyODEFunctionGenerator(right_hand_side,
                                                 sys.coordinates,
                                                 sys.speeds,
                                                 constants,
                                                 constants_arg_type=p_arg_type,
                                                 specifieds_arg_type=r_arg_type)

                assert g.specifieds_arg_type is None

                rhs = g.generate()

                xdot = rhs(x, 0.0, p[p_arg_type])

                try:
                    np.testing.assert_allclose(xdot, last_xdot)
                except NameError:
                    pass

                last_xdot = xdot
