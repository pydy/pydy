#!/usr/bin/env python

from sympy import symbols, trigsimp, zeros, Matrix
from sympy.physics.mechanics import *

"""
The disc spins about its Y axis and rolls on a plane which is normal to
the inertial reference frame's Z axis.

The essential generalized coordinates are:

q1 : yaw
q2 : roll
q3 : pitch

Two ignorable coordinates will be added to help create the visualization:

q4 : x ground contact point
q5 : y ground contant point

The generalized speeds are defined:

u1 : q2'
u2 : sin(q2) * q1'
u3 : cos(q2) * q1'

The model constants are:

r : radius of the disc
m : mass of the disc
g : acceleration due to gravity

"""

mechanics_printing(pretty_print=True)

q1, q2, q3, q4, q5 = dynamicsymbols('q1 q2 q3 q4 q5')
u1, u2, u3 = dynamicsymbols('u1 u2 u3')

q1d, q2d, q3d, q4d, q5d = dynamicsymbols('q1 q2 q3 q4 q5', 1)
u1d, u2d, u3d = dynamicsymbols('u1 u2 u3', 1)

r, m, g = symbols('r m g')

# The disc is rotated relative to the inertial reference frame by body fixed
# ZXY rotations.
N = ReferenceFrame('N')
Y = N.orientnew('Y', 'Axis', [q1, N.z])  # yaw angle
L = Y.orientnew('L', 'Axis', [q2, Y.x])  # lean/roll angle
R = L.orientnew('R', 'Axis', [q3, L.y])  # pitch angle

# The kinematic differential equations are defined such that the generalized
# speeds are the measure numbers of the lean frame's angular velocity.
w_R_N_qd = R.ang_vel_in(N)
R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)
kd = [dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L]

# C is the contact point and Ro is the mass center of the disc.
C = Point('C')
C.set_vel(N, 0)

Ro = C.locatenew('Ro', r * L.z)
Ro.v2pt_theory(C, N, R)

# The inertia of the disc can be defined in the L frame, instead of the R
# frame, because the inertia in the R frame is constant with respect to the
# L frame.
I = inertia(L,
            m / 4 * r**2,
            m / 2 * r**2,
            m / 4 * r**2)

# The force due to gravity is the only non-contributing force in the system
# and it acts in the negative inertial Z direction.
ForceList = [(Ro, -m * g * Y.z)]

# There is one rigid body in the system.
BodyR = RigidBody('BodyR', Ro, R, m, (I, Ro))
BodyList = [BodyR]

# Form the equations of motion.
kane = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3], kd_eqs=kd)
fr, frstar = kane.kanes_equations(ForceList, BodyList)

# Now we augment the kinematical equations of motion with the differential
# equations that govern the translation of the disc along it's rolling path
# so that we can integrate them to find the disc's path in inertial space.
# These will be used to make the visualization more interesting and
# realistic.

# This unit vector is always tagent to the path on the ground that the
# contact point traces out.
ground_path_tangent = cross(R.y, N.z).normalize()

# The magnitude of the ground contact velocity vector is defined by how fast
# the disc is rolling. The angular speed of the disc is defined as u2 and
# the ground contact speed is found by multiplying by the radius of the
# disc.
ground_path_velocity = u2 * r * ground_path_tangent

# A point that gives the translating center of mass is then:
translating_com = C.locatenew('Ro_t', q4 * N.x + q5 * N.y + r * L.z)

# The forcing vector can be simplified greatly:
F = Matrix(trigsimp(kane.forcing_full))

# Now the mass matrix and forcing vector can be augmented by the new
# kinematical differential equations. We simply insert the two new equations
# in after the existing kinematical equations.
M = Matrix(kane.mass_matrix_full)
M = M.col_insert(3, zeros(6, 2))
M = M.row_insert(3, Matrix([0, 0, 0, 1, 0, 0, 0, 0]).reshape(1, 8))
M = M.row_insert(4, Matrix([0, 0, 0, 0, 1, 0, 0, 0]).reshape(1, 8))

F = F.row_insert(3, Matrix([ground_path_velocity.dot(N.x)]))
F = F.row_insert(4, Matrix([ground_path_velocity.dot(N.y)]))
