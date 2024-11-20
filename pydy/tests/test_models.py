#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

# external libraries
import sympy as sm
import sympy.physics.mechanics as me
import numpy as np

# local
from ..models import multi_mass_spring_damper, n_link_pendulum_on_cart

def test_multi_mass_spring_damper():

    k0, m0, g, c0 = sm.symbols('k0, m0, g, c0')
    x0, v0, f0 = me.dynamicsymbols('x0, v0, f0')

    sys = multi_mass_spring_damper()

    assert sys.constants_symbols == {k0, c0, m0}
    assert sys.specifieds_symbols == set()
    assert sys.coordinates == [x0]
    assert sys.speeds == [v0]
    assert sys.states == [x0, v0]

    expected_mass_matrix_full = sm.Matrix([[1, 0], [0, m0]])
    expected_forcing_full = sm.Matrix([[v0],
                                       [-c0 * v0 - k0 * x0]])

    assert sm.simplify(sys.eom_method.mass_matrix_full -
                       expected_mass_matrix_full) == sm.zeros(2, 2)
    assert sm.simplify(sys.eom_method.forcing_full -
                       expected_forcing_full) == sm.zeros(2, 1)


def test_multi_mass_spring_damper_with_forces():

    k0, m0, g, c0 = sm.symbols('k0, m0, g, c0')
    x0, v0, f0 = me.dynamicsymbols('x0, v0, f0')

    sys = multi_mass_spring_damper(apply_gravity=True,
                                   apply_external_forces=True)

    assert sys.constants_symbols == {k0, m0, g, c0}
    assert sys.specifieds_symbols == {f0}
    assert sys.coordinates == [x0]
    assert sys.speeds == [v0]
    assert sys.states == [x0, v0]

    expected_mass_matrix_full = sm.Matrix([[1, 0], [0, m0]])
    expected_forcing_full = sm.Matrix([[v0],
                                       [-c0 * v0 + g * m0 - k0 * x0 + f0]])

    assert sm.simplify(sys.eom_method.mass_matrix_full -
                       expected_mass_matrix_full) == sm.zeros(2, 2)
    assert sm.simplify(sys.eom_method.forcing_full -
                       expected_forcing_full) == sm.zeros(2, 1)


def test_multi_mass_spring_damper_double():

    m0, m1, c0, c1, k0, k1, g = sm.symbols('m0, m1, c0, c1, k0, k1, g')
    x0, x1, v0, v1, f0, f1 = me.dynamicsymbols('x0, x1, v0, v1, f0, f1')

    sys = multi_mass_spring_damper(2, True, True)

    assert sys.constants_symbols == \
        {c1, m1, k0, c0, k1, g, m0}
    assert sys.specifieds_symbols == {f0, f1}
    assert sys.coordinates == [x0, x1]
    assert sys.speeds == [v0, v1]
    assert sys.states == [x0, x1, v0, v1]

    expected_mass_matrix_full = sm.Matrix([[1, 0,       0,  0],
                                           [0, 1,       0,  0],
                                           [0, 0, m0 + m1, m1],
                                           [0, 0,      m1, m1]])
    expected_forcing_full = sm.Matrix([[v0],
                                       [v1],
                                       [-c0 * v0 + g * m0 + g * m1 - k0 * x0 + f0 + f1],
                                       [-c1 * v1 + g * m1 - k1 * x1 + f1]])

    assert sm.simplify(sys.eom_method.mass_matrix_full -
                       expected_mass_matrix_full) == sm.zeros(4, 4)
    assert sm.simplify(sys.eom_method.forcing_full -
                       expected_forcing_full) == sm.zeros(4, 1)


def test_n_link_pendulum_on_cart_regression():
    num_time_steps = 1e7
    n = 3
    method = 'lambdify'

    parameter_vals = [9.81, 0.01 / n]
    m = sm.symbols('m:{}'.format(n + 1))
    l = sm.symbols('l:{}'.format(n))
    g = sm.symbols('g')

    sys = n_link_pendulum_on_cart(n, cart_force=False)
    sys.times = np.linspace(0, 3, num_time_steps)

    x0 = np.hstack(
        (0,
         np.pi/2*np.ones(len(sys.coordinates) - 1),
         1e-3*np.ones(len(sys.speeds))))
    sys.initial_conditions = dict(zip(sys.states, x0))

    constants = [g, m[0]]
    for i in range(n):
        constants += [l[i], m[i + 1]]

    sys.constants = dict(zip(constants, np.array(parameter_vals)))
    sys.generate_ode_function(generator=method)
    x = sys.integrate()
    #np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)),
    #                     'n_link_pendulum_on_cart_regression.npy'), x[::1000])

    datafile = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'n_link_pendulum_on_cart_regression.npy')
    expected_x = np.load(datafile)

    rtol = 1e-7
    atol = 1e-7
    x = x[::1000] # compare every 1000th sample
    np.testing.assert_allclose(x[:, :n+1], expected_x[:, :n+1], rtol, atol)

    # allow greater tolerance in speeds
    np.testing.assert_allclose(x[:, n+1:], expected_x[:, n+1:],
                               rtol*100, atol*100)
