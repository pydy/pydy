#This one is compatible with SymPy 0.7.4.1

from sympy import symbols, trigsimp, zeros, Matrix
from sympy.physics.mechanics import *

"""The disc spins about its Y axis and rolls on a plane which is normal to
the interial reference frame's Z axis."""

mechanics_printing()

"""
q1 : pitch
q2 : roll
q3 : yaw
q4 : x ground contact point
q5 : y ground contant point

r : radius of the disc
m : mass of the disc
g : acceleration due to gravity
"""

q1, q2, q3, q4, q5, u1, u2, u3 = dynamicsymbols('q1 q2 q3 q4 q5 u1 u2 u3')
q1d, q2d, q3d, q4d, q5d, u1d, u2d, u3d = dynamicsymbols('q1 q2 q3 q4 q5 u1 u2 u3', 1)
r, m, g = symbols('r m g')

# The disc is rotated relative to the inertial reference frame by body fixed
# YXZ rotations.
N = ReferenceFrame('N')
Y = N.orientnew('Y', 'Axis', [q1, N.y])  # pitch angle
L = Y.orientnew('L', 'Axis', [q2, Y.x])  # roll angle
R = L.orientnew('R', 'Axis', [q3, L.z])  # yaw angle

w_R_N_qd = R.ang_vel_in(N)
R.set_ang_vel(N, u1 * L.x + u2 * L.y + u3 * L.z)
kd = [dot(R.ang_vel_in(N) - w_R_N_qd, uv) for uv in L]

# C is the contact point and Dmc is the mass center of the disc.
C = Point('C')
C.set_vel(N, 0)
Dmc = C.locatenew('Dmc', r * L.z)
Dmc.v2pt_theory(C, N, R)
I = inertia(L, m / 4 * r**2, m / 2 * r**2, m / 4 * r**2)
mprint(I)

# The force due to gravity is the only non-contributing force in the system
# and it acts in the negative inertial Z direction.
ForceList = [(Dmc, -m * g * Y.z)]


BodyD = RigidBody('BodyD', Dmc, R, m, (I, Dmc))
BodyList = [BodyD]


kane = KanesMethod(N, q_ind=[q1, q2, q3], u_ind=[u1, u2, u3], kd_eqs=kd)
fr, frstar = kane.kanes_equations(ForceList, BodyList)

# Now we augment the kinematical equations of motion with the differential
# equations that govern the translation of the disc along it's rolling path.

# This unit vector is always tagent to the path on the ground that the
# contact point follows.
ground_path_tangent = cross(R.y, N.z).normalize()

# The magnitude of the ground contact velocity vector is defined by how fast
# the disc is rolling. The angular speed of the disc is found by
angular_speed = dot(R.ang_vel_in(N), R.y)

ground_path_velocity = angular_speed * r * ground_path_tangent

translating_contact = C.locatenew('C2', q4 * N.x + q5 * N.y)
translating_com = translating_contact.locatenew('D2', r * L.z)

kd = kane.kindiffdict()

kd[q4d] = ground_path_velocity.dot(N.x)
kd[q5d] = ground_path_velocity.dot(N.y)

F = Matrix(trigsimp(kane.forcing_full))

F = F.row_insert(3, Matrix([kd[q4d]]))
F = F.row_insert(4, Matrix([kd[q5d]]))

M = Matrix(kane.mass_matrix_full)
M = M.col_insert(3, zeros(6, 2))
M = M.row_insert(3, Matrix([0, 0, 0, 1, 0, 0, 0, 0]).reshape(1, 8))
M = M.row_insert(4, Matrix([0, 0, 0, 0, 1, 0, 0, 0]).reshape(1, 8))
