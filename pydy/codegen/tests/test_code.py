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
from .models import generate_mass_spring_damper_equations_of_motion


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

class TestCode():

    def test_generate_ode_function(self):

        m, k, c, g, F, x, v = np.random.random(7)

        args = {'constants': np.array([m, k, c, g]),
                'specified': np.array([F])}

        states = np.array([x, v])

        mass_matrix, forcing_vector, constants, coordinates, speeds, specified = \
            generate_mass_spring_damper_equations_of_motion()

        expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x + F)])

        backends = ['lambdify', 'theano', 'cython']

        for backend in backends:
            rhs = generate_ode_function(mass_matrix, forcing_vector,
                                        constants, coordinates, speeds,
                                        specified, generator=backend)
            dx = rhs(states, 0.0, args)

            testing.assert_allclose(dx, expected_dx)

        # Now try it with a function defining the specified quantities.
        args['specified'] = lambda x, t: np.sin(t)

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
