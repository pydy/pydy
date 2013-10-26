#!/usr/bin/env python

from sympy import symbols
import sympy.physics.mechanics as me
from code import numeric_right_hand_side, generate_mass_forcing_cython_code

"""
This is a simple mass, spring, damper system with gravity and an external
force.

 / / / / / / / /
----------------
   |    |     |   | g
   \   | |    |   V
 k /   --- c  |
   |    |     | x, v
  --------    V
  |  m   | -----
  --------
     | F
     V

"""

# derviation
mass, stiffness, damping, gravity = symbols('m, k, c, g')

position, speed = me.dynamicsymbols('x v')
positiond = me.dynamicsymbols('x', 1)
force = me.dynamicsymbols('F')

ceiling = me.ReferenceFrame('N')

origin = me.Point('origin')
origin.set_vel(ceiling, 0)

center = origin.locatenew('center', position * ceiling.x)
center.set_vel(ceiling, speed * ceiling.x)

block = me.Particle('block', center, mass)

kinematic_equations = [speed - positiond]
forces = [(center, (force + mass * gravity - stiffness * position - damping
                    * speed) * ceiling.x)]
bodies = [block]

kane = me.KanesMethod(ceiling, q_ind=[position], u_ind=[speed],
                      kd_eqs=kinematic_equations)
fr, frstar = kane.kanes_equations(forces, bodies)

# TODO: generate code for some output quantities too
kinetic_energy = block.kinetic_energy(ceiling)

# important parameters
coordinates = (position,)
speeds = (speed,)
specified = (force,)
constants = (mass, stiffness, damping, gravity)

# TODO : We may need to substitute in the u's before passing the results of
# Kane into the code generators.
kinematic_map = {positiond: speed}

def test_c_code_gen():
    prefix = 'actual_mass_forcing'
    generate_mass_forcing_cython_code(prefix, kane.mass_matrix_full,
                                      kane.forcing_full,
                                      coordinates,
                                      speeds,
                                      specified,
                                      constants)
    # TODO : Generate the shared library, compute the mass and forcing
    # matrices, and assert against expected.

parameter_values = [10.0, 5.0, 0.1, 9.8]
right_hand_side = numeric_right_hand_side(kane, constants, [kinetic_energy])
right_hand_side([1.0, 5.0], 0.1, parameter_values)
