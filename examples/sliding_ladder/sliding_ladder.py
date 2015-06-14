#!/usr/bin/env python

"""Sliding Ladder in a Frictionless World

A ladder is placed on a floor and is leaning against a wall. Both of the
surfaces in contact with the ladder are frictionless. There are two reaction
forces on the ladder (one from the floor and one from the wall). The ladder
falls under the influence of gravity.

We need only one angle to fully define the position of the ladder. The mass,
length, and moment of intertia of the ladder are given.

"""

import sympy as sm
import sympy.physics.mechanics as me

# First we define all the symbols we need: the angle, q, and angular
# velocity, u, of the ladder.
q, u = me.dynamicsymbols('q u')

# The time derivatives of the symbols declared above:
qd, ud = me.dynamicsymbols('q u', 1)

# The mass, length, and moment of inertia of the ladder and the acceleration
# due to gravity.
m, l, g, Izz = sm.symbols('m l g Izz')

# We define the inertial frame and a lean frame for the ladder.
N = me.ReferenceFrame('N')
L = N.orientnew('L', 'Axis', [q, N.z])

# And then set the angular velocity of the lean frame in the inertial frame.
# The angular acceleration will be automatically computed when needed.
L.set_ang_vel(N, u * N.z)

# Now the origin and the center of mass of the ladder are defined.
O = me.Point('O')
A = me.Point('A')

# And we use the length and angle to locate the center of mass relative to
# the origin.
A.set_pos(O, -l / 2 * sm.cos(q) * N.x + l / 2 * sm.sin(q) * N.y)

# Take the derivative of the position of the center of mass to get the
# corresponding velocity.
O.set_vel(N, 0)
A.set_vel(N, l / 2 * u * sm.sin(q) * N.x + l / 2 * u * sm.cos(q) * N.y)

# The ladder can now be defined as a rigid body.
ladder = me.RigidBody('ladder', A, L, m, (me.inertia(L, 0, 0, Izz), A))

# Set up all the inputs to Kane's method.
kd = [u - qd]
body_list = [ladder]
force_list = [(A, -m * g * N.y)]

# Finally we solve the dynamics
kane = me.KanesMethod(N, q_ind=[q], u_ind=[u], kd_eqs=kd)
kane.kanes_equations(force_list, body_list)

# The mass matrix and the forcing function can be taken out of the
# KanesMethod object.
MM = kane.mass_matrix
forcing = kane.forcing

# And those can be used to find the equations of motion.
rhs = MM.inv() * forcing
kdd = kane.kindiffdict()

rhs = rhs.subs(kdd)
rhs.simplify()
print(rhs)
