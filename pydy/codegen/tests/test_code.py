#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
import sympy as sm
import numpy as np
from numpy import testing

# local libraries
from ..code import generate_ode_function
from ...models import multi_mass_spring_damper, n_link_pendulum_on_cart

# TODO : Remove these tests before PyDy 0.4.0. They are just here to make
# sure old stuff still runs.


class TestCodeRHSArgs():

    def test_rhs_args(self):

        sys = n_link_pendulum_on_cart(3, True, True)

        constants = sys.constants_symbols
        specifieds = tuple(sys.specifieds_symbols)

        args = (sys.eom_method.mass_matrix_full,
                sys.eom_method.forcing_full,
                constants,
                sys.coordinates,
                sys.speeds)

        kwargs = {'specified': specifieds}

        rhs = generate_ode_function(*args, **kwargs)

        x = np.array(np.random.random(len(sys.coordinates + sys.speeds)))

        rhs_args = {'constants': {k: 1.0 for k in constants}}

        rhs_args['specified'] = dict(zip(specifieds, [1.0, 2.0, 3.0, 4.0]))

        xd_01 = rhs(x, 0.0, rhs_args)

        rhs_args['specified'] = {
            specifieds: lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])}

        xd_02 = rhs(x, 0.0, rhs_args)

        # There are four specified inputs available.
        rhs_args['specified'] = {
            specifieds[0]: lambda x, t: np.ones(1),
            (specifieds[3], specifieds[1]): lambda x, t: np.array([4.0, 2.0]),
            specifieds[2]: 3.0 * np.ones(1)}

        xd_03 = rhs(x, 0.0, rhs_args)

        testing.assert_allclose(xd_01, xd_02)
        testing.assert_allclose(xd_01, xd_03)

        # Test old and efficient RHS args.
        rhs_args['specified'] = np.array([1.0, 2.0, 3.0, 4.0])
        xd_04 = rhs(x, 0.0, rhs_args)
        testing.assert_allclose(xd_01, xd_04)

        rhs_args['specified'] = lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])
        xd_05 = rhs(x, 0.0, rhs_args)
        testing.assert_allclose(xd_01, xd_05)


class TestCode():

    def test_generate_ode_function(self):

        system = multi_mass_spring_damper(1, True, True)

        args = (system.eom_method.mass_matrix_full,
                system.eom_method.forcing_full,
                system.constants_symbols,
                system.coordinates,
                system.speeds)

        kwargs = {'specified': tuple(system.specifieds_symbols)}

        mass_matrix, forcing_vector, constants, coordinates, speeds = args
        specifieds = kwargs['specified']

        F, x, v = np.random.random(3)

        states = np.array([x, v])

        constants_map = dict(zip(constants, np.random.random(len(constants))))

        m = constants_map[sm.symbols('m0')]
        k = constants_map[sm.symbols('k0')]
        c = constants_map[sm.symbols('c0')]
        g = constants_map[sm.symbols('g')]

        expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x + F)])

        rhs_args = {'constants': constants_map,
                    'specified': {specifieds[0]: F}}

        backends = ['lambdify', 'theano', 'cython']

        for backend in backends:
            rhs = generate_ode_function(mass_matrix, forcing_vector,
                                        constants, coordinates, speeds,
                                        specifieds, generator=backend)
            dx = rhs(states, 0.0, rhs_args)

            testing.assert_allclose(dx, expected_dx)

        # Now try it with a function defining the specified quantities.
        rhs_args['specified'] = {specifieds[0]: lambda x, t: np.sin(t)}

        t = 14.345

        expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x +
                                              np.sin(t))])

        for backend in backends:
            rhs = generate_ode_function(mass_matrix, forcing_vector,
                                        constants, coordinates, speeds,
                                        specifieds, generator=backend)
            dx = rhs(states, t, rhs_args)

            testing.assert_allclose(dx, expected_dx)

        # Now try it without specified values.
        system = multi_mass_spring_damper(1, True)

        args = (system.eom_method.mass_matrix_full,
                system.eom_method.forcing_full,
                system.constants_symbols,
                system.coordinates,
                system.speeds)

        mass_matrix, forcing_vector, constants, coordinates, speeds = args

        expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x)])

        for backend in backends:
            rhs = generate_ode_function(mass_matrix, forcing_vector,
                                        constants, coordinates, speeds,
                                        specifieds, generator=backend)
            dx = rhs(states, 0.0, rhs_args)

            testing.assert_allclose(dx, expected_dx)
