#!/usr/bin/env python

# external libraries
from code import numeric_right_hand_side, generate_mass_forcing_cython_code

# local libraries
from models import generate_mass_spring_damper_equations_of_motion


def test_cython_code_gen():

    mass_matrix, forcing_vector, constants, coordinates, speeds, specified = \
        generate_mass_spring_damper_equations_of_motion()

    prefix = 'actual_mass_forcing'
    generate_mass_forcing_cython_code(prefix, mass_matrix, forcing_vector,
                                      constants, coordinates, speeds,
                                      specified)

    # TODO : Generate the shared library, compute the mass and forcing
    # matrices, and assert against expected.


def test_numeric_evaluation():

    parameter_values = [10.0, 5.0, 0.1, 9.8]
    right_hand_side = numeric_right_hand_side(kane, constants)
    right_hand_side([1.0, 5.0], 0.1, parameter_values)
