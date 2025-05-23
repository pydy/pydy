#!/usr/bin/env python

from random import choice
import warnings

import numpy as np
import scipy as sp
import sympy as sm
from pydy.codegen.ode_function_generators import generate_ode_function

Cython = sm.external.import_module('Cython')
theano = sm.external.import_module('theano')

from ... import models
from ..ode_function_generators import (ODEFunctionGenerator,
                                       LambdifyODEFunctionGenerator,
                                       CythonODEFunctionGenerator,
                                       TheanoODEFunctionGenerator)

from ...utils import PyDyImportWarning

warnings.simplefilter('once', PyDyImportWarning)


def test_symbolic_lusolve_full_mass_matrix():
    sys = models.n_link_pendulum_on_cart(n=5, cart_force=False,
                                         joint_torques=False)

    g_symbolic_solve = CythonODEFunctionGenerator(
        sys.eom_method.forcing_full,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix_full,
        linear_sys_solver='sympy')
    rhs_symbolic_solve = g_symbolic_solve.generate()

    g_numeric_solve = CythonODEFunctionGenerator(
        sys.eom_method.forcing_full,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix_full,
        linear_sys_solver='numpy')
    rhs_numeric_solve = g_numeric_solve.generate()

    x = np.random.random(g_symbolic_solve.num_coordinates +
                         g_symbolic_solve.num_speeds)
    t = 5.125
    p = np.random.random(g_symbolic_solve.num_constants)

    np.testing.assert_allclose(rhs_numeric_solve(x, t, p),
                               rhs_symbolic_solve(x, t, p))


def test_symbolic_lusolve_min_mass_matrix():
    sys = models.n_link_pendulum_on_cart(n=5, cart_force=False,
                                         joint_torques=False)
    kin_diff_eqs = sys.eom_method.kindiffdict()
    coord_derivs = sm.Matrix([kin_diff_eqs[c.diff()] for c in
                              sys.coordinates])

    g_symbolic_solve = CythonODEFunctionGenerator(
        sys.eom_method.forcing,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix,
        coordinate_derivatives=coord_derivs,
        linear_sys_solver='sympy')
    rhs_symbolic_solve = g_symbolic_solve.generate()

    g_numeric_solve = CythonODEFunctionGenerator(
        sys.eom_method.forcing,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix,
        coordinate_derivatives=coord_derivs,
        linear_sys_solver='numpy')
    rhs_numeric_solve = g_numeric_solve.generate()

    x = np.random.random(g_symbolic_solve.num_coordinates +
                         g_symbolic_solve.num_speeds)
    t = 5.125
    p = np.random.random(g_symbolic_solve.num_constants)

    np.testing.assert_allclose(rhs_numeric_solve(x, t, p),
                               rhs_symbolic_solve(x, t, p))


def test_cse_same_numerical_results():
    # NOTE : This ensurses that the same results are always given for the sympy
    # cse outputs, which seem to change every version.

    if not Cython:
        return

    sys = models.n_link_pendulum_on_cart(n=5, cart_force=False,
                                         joint_torques=False)

    g_no_cse = CythonODEFunctionGenerator(
        sys.eom_method.forcing_full,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix_full,
        cse=False)
    rhs_func_no_cse = g_no_cse.generate()

    g_cse = CythonODEFunctionGenerator(
        sys.eom_method.forcing_full,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix_full,
        cse=True)
    rhs_func_cse = g_cse.generate()

    g_lam_cse = LambdifyODEFunctionGenerator(
        sys.eom_method.forcing_full,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix_full,
        cse=True)
    rhs_func_lam_cse = g_lam_cse.generate()

    g_lam_no_cse = LambdifyODEFunctionGenerator(
        sys.eom_method.forcing_full,
        sys.coordinates,
        sys.speeds,
        sys.constants_symbols,
        mass_matrix=sys.eom_method.mass_matrix_full,
        cse=False)
    rhs_func_lam_no_cse = g_lam_no_cse.generate()

    x = np.random.random(g_cse.num_coordinates + g_cse.num_speeds)
    t = 5.125
    p = np.random.random(g_cse.num_constants)

    np.testing.assert_allclose(rhs_func_no_cse(x, t, p),
                               rhs_func_cse(x, t, p))

    np.testing.assert_allclose(rhs_func_lam_no_cse(x, t, p),
                               rhs_func_lam_cse(x, t, p))

    np.testing.assert_allclose(rhs_func_lam_cse(x, t, p),
                               rhs_func_cse(x, t, p))


class TestODEFunctionGenerator(object):

    def setup_method(self):

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

        assert g.linear_sys_solver == 'numpy'
        assert g._solve_linear_system == np.linalg.solve

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver='numpy')

        assert g.linear_sys_solver == 'numpy'
        assert g._solve_linear_system == np.linalg.solve

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver='scipy')

        assert g.linear_sys_solver == 'scipy'
        assert g._solve_linear_system == sp.linalg.solve

        solver = lambda A, b: np.dot(np.inv(A), b)

        g = ODEFunctionGenerator(self.sys.eom_method.forcing_full,
                                 self.sys.coordinates,
                                 self.sys.speeds,
                                 self.sys.constants_symbols,
                                 mass_matrix=self.sys.eom_method.mass_matrix_full,
                                 linear_sys_solver=solver)

        assert g.linear_sys_solver == solver
        assert g._solve_linear_system == solver

    def test_no_constants(self):
        sys = models.multi_mass_spring_damper()
        constant_vals = {sm.Symbol('m0'): 1.0, sm.Symbol('c0'): 2.0, sm.Symbol('k0'): 3.0}

        sym_rhs = sys.eom_method.rhs()

        rhs = generate_ode_function(sym_rhs, sys.coordinates, sys.speeds, constant_vals)
        rhs2 = generate_ode_function(sym_rhs.subs(constant_vals), sys.coordinates, sys.speeds)

        assert np.array_equal(rhs(np.array([1.0, 2.0]), 0.0, constant_vals),
                              rhs2(np.array([1.0, 2.0]), 0.0))


class TestODEFunctionGeneratorSubclasses(object):

    ode_function_subclasses = [LambdifyODEFunctionGenerator]

    if Cython:
        ode_function_subclasses.append(CythonODEFunctionGenerator)
    else:
        warnings.warn("Cython was not found so the related tests are being"
                      " skipped.", PyDyImportWarning)

    if theano:
        ode_function_subclasses.append(TheanoODEFunctionGenerator)
    else:
        warnings.warn("Theano was not found so the related tests are being"
                      " skipped.", PyDyImportWarning)

    def setup_method(self):

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

        x = np.random.random(len(sys.states))
        t_val = 1.23

        r_array = t_val*x[:4]
        r_dct_1 = dict(zip(specifieds, r_array))
        r_dct_2 = {tuple(specifieds): lambda x, t: t*x[:4]}
        r_dct_3 = {specifieds[0]: lambda x, t: t*x[0],
                   (specifieds[3], specifieds[1]): lambda x, t: t*x[[3, 1]],
                   specifieds[2]: r_array[2]}
        r_func = lambda x, t: t*x[:4]

        r = {}
        r[None] = choice([r_array, r_dct_1, r_dct_2, r_dct_3, r_func])
        r['array'] = r_array
        r['dictionary'] = choice([r_dct_1, r_dct_2, r_dct_3])
        r['function'] = r_func

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

                xdot = rhs(x, t_val, r[r_arg_type], p[p_arg_type])

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

                xdot = rhs(x, t_val, p[p_arg_type])

                try:
                    np.testing.assert_allclose(xdot, last_xdot)
                except NameError:
                    pass

                last_xdot = xdot

    def test_rhs_docstring(self):

        sys = models.n_link_pendulum_on_cart(2, False, False)
        right_hand_side = sys.eom_method.rhs()

        constants = list(sm.ordered(sys.constants_symbols))

        constants_arg_types = [None, 'array', 'dictionary']

        rhs_doc_template = \
"""\
Returns the derivatives of the states, i.e. numerically evaluates the right
hand side of the first order differential equation.

x' = f(x, t,{specified_call_sig} p)

Parameters
==========
x : ndarray, shape(6,)
    The state vector is ordered as such:
        - q0(t)
        - q1(t)
        - q2(t)
        - u0(t)
        - u1(t)
        - u2(t)
t : float
    The current time.{specifieds_explanation}{constants_explanation}
Returns
=======
dx : ndarray, shape(6,)
    The derivative of the state vector.

"""

        constants_doc_templates = {}

        constants_doc_templates['dictionary'] = \
"""
p : dictionary len(6)
    A dictionary that maps the constants symbols to their numerical values
    with at least these keys:
        - g
        - l0
        - l1
        - m0
        - m1
        - m2
"""

        constants_doc_templates['array'] = \
"""
p : ndarray shape(6,)
    A ndarray of floats that give the numerical values of the constants in
    this order:
        - g
        - l0
        - l1
        - m0
        - m1
        - m2
"""

        constants_doc_templates[None] = \
"""
p : dictionary len(6) or ndarray shape(6,)
    Either a dictionary that maps the constants symbols to their numerical
    values or an array with the constants in the following order:
        - g
        - l0
        - l1
        - m0
        - m1
        - m2
"""

        for p_arg_type in constants_arg_types:

            _rhs_doc_template = rhs_doc_template.format(**{
                'specified_call_sig': '',
                'specifieds_explanation': '',
                'constants_explanation': constants_doc_templates[p_arg_type]
                })

            g = LambdifyODEFunctionGenerator(right_hand_side,
                                             sys.coordinates,
                                             sys.speeds,
                                             constants,
                                             constants_arg_type=p_arg_type)

            rhs = g.generate()

            assert (_rhs_doc_template == rhs.__doc__)

        sys = models.n_link_pendulum_on_cart(2, True, True)
        right_hand_side = sys.eom_method.rhs()

        constants = list(sm.ordered(sys.constants_symbols))
        specifieds = list(sm.ordered(sys.specifieds_symbols))

        specifieds_arg_types = [None, 'array', 'function', 'dictionary']

        specifieds_doc_templates = {}

        specifieds_doc_templates[None] = \
"""
r : dictionary; ndarray, shape(3,); function

    There are three options for this argument. (1) is more flexible but
    (2) and (3) are much more efficient.

    (1) A dictionary that maps the specified functions of time to floats,
    ndarrays, or functions that produce ndarrays. The keys can be a single
    specified symbolic function of time or a tuple of symbols. The total
    number of symbols must be equal to 3. If the value is a
    function it must be of the form g(x, t), where x is the current state
    vector ndarray and t is the current time float and it must return an
    ndarray of the correct shape. For example::

      r = {a: 1.0,
           (d, b) : np.array([1.0, 2.0]),
           (e, f) : lambda x, t: np.array(x[0], x[1]),
           c: lambda x, t: np.array(x[2])}

    (2) A ndarray with the specified values in the correct order and of the
    correct shape.

    (3) A function that must be of the form g(x, t), where x is the current
    state vector and t is the current time and it must return an ndarray of
    the correct shape.

    The specified inputs are, in order:
        - F(t)
        - T1(t)
        - T2(t)"""

        specifieds_doc_templates['array'] = \
"""
r : ndarray, shape(3,)

    A ndarray with the specified values in the correct order and of the
    correct shape.

    The specified inputs are, in order:
        - F(t)
        - T1(t)
        - T2(t)"""

        specifieds_doc_templates['dictionary'] = \
"""
r : dictionary

    A dictionary that maps the specified functions of time to floats,
    ndarrays, or functions that produce ndarrays. The keys can be a single
    specified symbolic function of time or a tuple of symbols. The total
    number of symbols must be equal to 3. If the value is a
    function it must be of the form g(x, t), where x is the current state
    vector ndarray and t is the current time float and it must return an
    ndarray of the correct shape. For example::

      r = {a: 1.0,
           (d, b) : np.array([1.0, 2.0]),
           (e, f) : lambda x, t: np.array(x[0], x[1]),
           c: lambda x, t: np.array(x[2])}

    The specified inputs are, in order:
        - F(t)
        - T1(t)
        - T2(t)"""

        specifieds_doc_templates['function'] = \
"""
r : function

    A function that must be of the form g(x, t), where x is the current
    state vector and t is the current time and it must return an ndarray of
    shape(3,).

    The specified inputs are, in order:
        - F(t)
        - T1(t)
        - T2(t)"""

        for p_arg_type in constants_arg_types:
            for r_arg_type in specifieds_arg_types:

                _rhs_doc_template = rhs_doc_template.format(**{
                    'specified_call_sig': ' r,',
                    'specifieds_explanation': specifieds_doc_templates[r_arg_type],
                    'constants_explanation': constants_doc_templates[p_arg_type]
                    })

                g = LambdifyODEFunctionGenerator(right_hand_side,
                                                 sys.coordinates,
                                                 sys.speeds,
                                                 constants,
                                                 specifieds=specifieds,
                                                 constants_arg_type=p_arg_type,
                                                 specifieds_arg_type=r_arg_type)

                rhs = g.generate()

                assert (_rhs_doc_template == rhs.__doc__)
