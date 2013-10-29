#!/usr/bin/env python

# standard libary
import importlib
import os
import shutil
import glob

# external libraries
import numpy as np
from numpy import testing

# local libraries
from pydy_code_gen.code import (numeric_right_hand_side,
                                generate_mass_forcing_cython_code)
from models import generate_mass_spring_damper_equations_of_motion


def test_numeric_right_hand_side():

    m, k, c, g, F, x, v = np.random.random(7)

    args = {'constants': np.array([m, k, c, g]),
            'specified': np.array([F]),
            'num_coordinates': 1}

    states = np.array([x, v])

    mass_matrix, forcing_vector, constants, coordinates, speeds, specified = \
        generate_mass_spring_damper_equations_of_motion()

    expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x + F)])

    backends = ['lambdify', 'cython'] #, 'theano']
    for backend in backends:
        rhs = numeric_right_hand_side(mass_matrix, forcing_vector,
                                      constants, coordinates, speeds,
                                      specified, generator=backend)
        dx = rhs(states, 0.0, args)

        testing.assert_allclose(dx, expected_dx)

        if backend == 'cython':
            filelist = glob.glob('multibody_system*')
            for f in filelist:
                os.remove(f)
            shutil.rmtree('build')

    mass_matrix, forcing_vector, constants, coordinates, speeds, specified = \
        generate_mass_spring_damper_equations_of_motion(external_force=False)

    expected_dx = np.array([v, 1.0 / m * (-c * v + m * g - k * x)])

    # TODO : There is an import issue here, where the previous cython code
    # is loaded from the module and the import needs to be reloaded. Not
    # sure how to fix that.

    for backend in backends:
        rhs = numeric_right_hand_side(mass_matrix, forcing_vector,
                                      constants, coordinates, speeds,
                                      specified, generator=backend)
        dx = rhs(states, 0.0, args)

        testing.assert_allclose(dx, expected_dx)

        if backend == 'cython':
            filelist = glob.glob('multibody_system*')
            for f in filelist:
                os.remove(f)
            shutil.rmtree('build')


def test_generate_mass_forcing_cython_code():

    mass_matrix, forcing_vector, constants, coordinates, speeds, specified = \
        generate_mass_spring_damper_equations_of_motion()

    prefix = 'actual_mass_forcing'
    generate_mass_forcing_cython_code(prefix, mass_matrix, forcing_vector,
                                      constants, coordinates, speeds,
                                      specified)

    cython_module = importlib.import_module(prefix)

    # TODO : Generate the shared library, compute the mass and forcing
    # matrices, and assert against expected.
