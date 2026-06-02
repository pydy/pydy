#!/usr/bin/env python

from sympy import symbols
import sympy.physics.mechanics as me

print("Defining the problem.")

# The conical pendulum will have three links and three bobs.
n = 3

# Each link's orientation is described by two spaced fixed angles: alpha and
# beta.

# Generalized coordinates
alpha = me.dynamicsymbols('alpha:{}'.format(n))
beta = me.dynamicsymbols('beta:{}'.format(n))

# Generalized speeds
omega = me.dynamicsymbols('omega:{}'.format(n))
delta = me.dynamicsymbols('delta:{}'.format(n))

# At each joint there are point masses (i.e. the bobs).
m_bob = symbols('m:{}'.format(n))

# Each link is modeled as a cylinder so it will have a length, mass, and a
# symmetric inertia tensor.
l = symbols('l:{}'.format(n))
m_link = symbols('M:{}'.format(n))
Ixx = symbols('Ixx:{}'.format(n))
Iyy = symbols('Iyy:{}'.format(n))
Izz = symbols('Izz:{}'.format(n))

# Acceleration due to gravity will be used when prescribing the forces
# acting on the links and bobs.
g = symbols('g')

# Now defining an inertial reference frame for the system to live in. The Y
# axis of the frame will be aligned with, but opposite to, the gravity
# vector.

I = me.ReferenceFrame('I')

# Three reference frames will track the orientation of the three links.

A = me.ReferenceFrame('A')
A.orient(I, 'Space', [alpha[0], beta[0], 0], 'ZXY')

B = me.ReferenceFrame('B')
B.orient(A, 'Space', [alpha[1], beta[1], 0], 'ZXY')

C = me.ReferenceFrame('C')
C.orient(B, 'Space', [alpha[2], beta[2], 0], 'ZXY')

# Define the kinematical differential equations such that the generalized
# speeds equal the time derivative of the generalized coordinates.
kinematic_differentials = []
for i in range(n):
    kinematic_differentials.append(omega[i] - alpha[i].diff())
    kinematic_differentials.append(delta[i] - beta[i].diff())

# The angular velocities of the three frames can then be set.
A.set_ang_vel(I, omega[0] * I.z + delta[0] * I.x)
B.set_ang_vel(I, omega[1] * I.z + delta[1] * I.x)
C.set_ang_vel(I, omega[2] * I.z + delta[2] * I.x)

# The base of the pendulum will be located at a point O which is stationary
# in the inertial reference frame.
O = me.Point('O')
O.set_vel(I, 0)

# The location of the bobs (at the joints between the links) are created by
# specifiying the vectors between the points.
P1 = O.locatenew('P1', -l[0] * A.y)
P2 = P1.locatenew('P2', -l[1] * B.y)
P3 = P2.locatenew('P3', -l[2] * C.y)

# The velocities of the points can be computed by taking advantage that
# pairs of points are fixed on the referene frames.
P1.v2pt_theory(O, I, A)
P2.v2pt_theory(P1, I, B)
P3.v2pt_theory(P2, I, C)
points = [P1, P2, P3]

# Now create a particle to represent each bob.
Pa1 = me.Particle('Pa1', points[0], m_bob[0])
Pa2 = me.Particle('Pa2', points[1], m_bob[1])
Pa3 = me.Particle('Pa3', points[2], m_bob[2])
particles = [Pa1, Pa2, Pa3]

# The mass centers of each link need to be specified and, assuming a
# constant density cylinder, it is equidistance from each joint.
P_link1 = O.locatenew('P_link1', -l[0] / 2 * A.y)
P_link2 = P1.locatenew('P_link2', -l[1] / 2 * B.y)
P_link3 = P2.locatenew('P_link3', -l[2] / 2 * C.y)

# The linear velocities can be specified the same way as the bob points.
P_link1.v2pt_theory(O, I, A)
P_link2.v2pt_theory(P1, I, B)
P_link3.v2pt_theory(P2, I, C)

points_rigid_body = [P_link1, P_link2, P_link3]

# The inertia tensors for the links are defined with respect to the mass
# center of the link and the link's reference frame.
inertia_link1 = (me.inertia(A, Ixx[0], Iyy[0], Izz[0]), P_link1)
inertia_link2 = (me.inertia(B, Ixx[1], Iyy[1], Izz[1]), P_link2)
inertia_link3 = (me.inertia(C, Ixx[2], Iyy[2], Izz[2]), P_link3)

# Now rigid bodies can be created for each link.
link1 = me.RigidBody('link1', P_link1, A, m_link[0], inertia_link1)
link2 = me.RigidBody('link2', P_link2, B, m_link[1], inertia_link2)
link3 = me.RigidBody('link3', P_link3, C, m_link[2], inertia_link3)
links = [link1, link2, link3]

# The only contributing forces to the system is the force due to gravity
# acting on each particle and body.
forces = []

for particle in particles:
    mass = particle.mass
    point = particle.point
    forces.append((point, -mass * g * I.y))

for link in links:
    mass = link.mass
    point = link.masscenter
    forces.append((point, -mass * g * I.y))

# Make a list of all the particles and bodies in the system.
total_system = links + particles

# Lists of all generalized coordinates and speeds.
q = alpha + beta
u = omega + delta

# Now the equations of motion of the system can be formed.
print("Generating equations of motion.")
kane = me.KanesMethod(I, q_ind=q, u_ind=u, kd_eqs=kinematic_differentials)
fr, frstar = kane.kanes_equations(loads=forces, bodies=total_system)
print("Derivation complete.")
