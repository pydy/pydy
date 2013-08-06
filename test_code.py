#!/usr/bin/env python

from sympy import symbols
import sympy.physics.mechanics as me
from sympy.physics.mechanics.code import numeric_right_hand_side

mass, stiffness, damping, gravity = symbols('mass stiffness damping gravity')

y, v = me.dynamicsymbols('y v')
yd = me.dynamicsymbols('y', 1)

ground = me.ReferenceFrame('ground')
block_frame = ground.orientnew('B', 'Axis', [0, ground.z])

origin = me.Point('origin')
origin.set_vel(ground, 0)

center = origin.locatenew('center', y * ground.y)
center.set_vel(ground, v * ground.y)

block = me.Particle('block', center, mass)

kinematic_equations = [v - yd]
forces = [(center, mass * gravity * ground.y - stiffness * y * ground.y -
    damping * v * ground.y)]
bodies = [block]

kane = me.KanesMethod(ground, q_ind=[y], u_ind=[v], kd_eqs=kinematic_equations)
fr, frstar = kane.kanes_equations(forces, bodies)

kinetic_energy = block.kinetic_energy(ground)

parameters = [mass, stiffness, damping, gravity]
parameter_values = [10.0, 5.0, 0.1, 9.8]

right_hand_side = numeric_right_hand_side(kane, parameters, [kinetic_energy])

right_hand_side([1.0, 5.0], 0.1, parameter_values)
