#!/usr/bin/env python

# standard libary
import os
import shutil
import glob
import filecmp

# external libraries
import numpy as np
from numpy import testing

# local libraries
from ..code import generate_ode_function, CythonGenerator
from .models import (generate_mass_spring_damper_equations_of_motion,
                     generate_n_link_pendulum_on_cart_equations_of_motion)


class TestCythonGenerator():

    prefix = 'desired_mass_forcing'

    def test_write_cython_code(self):

        results = generate_mass_spring_damper_equations_of_motion()

        generator = CythonGenerator(self.prefix, *results)
        generator._write_cython_code()

        file_dir = os.path.split(__file__)[0]

        expected_files = ['desired_mass_forcing_c.c',
                          'desired_mass_forcing_c.h',
                          'desired_mass_forcing.pyx',
                          'desired_mass_forcing_setup.py']

        endings = ['_c.c', '_c.h', '.pyx', '_setup.py']

        for ending, expected_file in zip(endings, expected_files):
            created = self.prefix + ending
            expected = os.path.join(file_dir, 'expected_cython',
                                    expected_file)
            assert filecmp.cmp(created, expected)

    def teardown(self):

        # clean up the cython crud
        files = glob.glob(self.prefix + '*')
        for f in files:
            os.remove(f)


class TestCodeRHSArgs():

    def test_rhs_args(self):

        sys = generate_n_link_pendulum_on_cart_equations_of_motion(3,
                                                                   True,
                                                                   True)
        constants = sys[2]
        coordinates = sys[3]
        speeds = sys[4]
        specified = sys[5]

        rhs = generate_ode_function(*sys)

        x = np.array(np.random.random(len(coordinates + speeds)))

        args = {'constants': {k: 1.0 for k in constants}}

        args['specified'] = dict(zip(specified, [1.0, 2.0, 3.0, 4.0]))

        xd_01 = rhs(x, 0.0, args)

        args['specified'] = {
                tuple(specified): lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])}

        xd_02 = rhs(x, 0.0, args)

        # There are four specified inputs available.
        args['specified'] = {specified[0]: lambda x, t: np.ones(1),
                (specified[3], specified[1]): lambda x, t: np.array([4.0, 2.0]),
                specified[2]: 3.0 * np.ones(1)}

        xd_03 = rhs(x, 0.0, args)

        testing.assert_allclose(xd_01, xd_02)
        testing.assert_allclose(xd_01, xd_03)

        # Test old and efficient RHS args.
        args['specified'] = np.array([1.0, 2.0, 3.0, 4.0])
        xd_04 = rhs(x, 0.0, args)
        testing.assert_allclose(xd_01, xd_04)

        args['specified'] = lambda x, t: np.array([1.0, 2.0, 3.0, 4.0])
        xd_05 = rhs(x, 0.0, args)
        testing.assert_allclose(xd_01, xd_05)

class TestCode():

    def test_generate_ode_function(self):

        system = generate_mass_spring_damper_equations_of_motion()
        mass_matrix = system[0]
        forcing_vector = system[1]
        constants = system[2]
        coordinates = system[3]
        speeds = system[4]
        specified = system[5]

        m, k, c, g, F, x, v = np.random.random(7)

        args = {'constants': dict(zip(constants, [m, k, c, g])),
                'specified': {specified[0]: F}}

        states = np.array([x, v])

        expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x + F)])

        backends = ['lambdify', 'theano', 'cython']

        for backend in backends:
            rhs = generate_ode_function(mass_matrix, forcing_vector,
                                        constants, coordinates, speeds,
                                        specified, generator=backend)
            dx = rhs(states, 0.0, args)

            testing.assert_allclose(dx, expected_dx)

        # Now try it with a function defining the specified quantities.
        args['specified'] = {specified[0]: lambda x, t: np.sin(t)}

        t = 14.345

        expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x +
                                              np.sin(t))])

        for backend in backends:
            rhs = generate_ode_function(mass_matrix, forcing_vector,
                                        constants, coordinates, speeds,
                                        specified, generator=backend)
            dx = rhs(states, t, args)

            testing.assert_allclose(dx, expected_dx)

        # Now try it without specified values.
        mass_matrix, forcing_vector, constants, coordinates, speeds, specified = \
            generate_mass_spring_damper_equations_of_motion(external_force=False)

        expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x)])

        for backend in backends:
            rhs = generate_ode_function(mass_matrix, forcing_vector,
                                        constants, coordinates, speeds,
                                        specified, generator=backend)
            dx = rhs(states, 0.0, args)

            testing.assert_allclose(dx, expected_dx)

    def teardown(self):

        # clean up the cython crud
        files = glob.glob('multibody_system*')
        for f in files:
            os.remove(f)
        shutil.rmtree('build')
