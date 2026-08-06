#!/usr/bin/env python

"""This example is setup to demonstrate several basic joint types on a
problem with no kinematic loops or nonholonomic constraints."""

import numpy as np
import sympy as sm
import sympy.physics.mechanics as me
from pydy.system import System
from pydy.viz import Plane, Cube, Sphere, VisualizationFrame, Scene

q = me.dynamicsymbols('q:8')
u = me.dynamicsymbols('u:8')

conical_length = sm.symbols('l_c')

plate_mass = sm.symbols('m_{pl}')
plate_length = sm.symbols('l_{pl}')

pendulum_length = sm.symbols('l_p')
pendulum_base_mass = sm.symbols('m_{p1}')
pendulum_bob_mass = sm.symbols('m_{p2}')

block_length = sm.symbols('l_b')
block_mass = sm.symbols('m_b')

k2, k3, k4, k5, k6 = sm.symbols('k2, k3, k4, k5, k6')

g = sm.symbols('g')

# Kinematical Differential Equations
# Make the simplest possible definitions for the kinematical differential
# equations: q' = u.

kin_diff_eqs = {q_i.diff(): u_i for q_i, u_i in zip(q, u)}

# Define the reference frames and their orientations.
# The conical pendulum uses quaternions as the generalized coordinates to
# avoid any singularities which would limit the rotation. The plate rotates
# around the conical pendulum shaft and the simple pendulum swings in the
# plate's plane.

newtonian_fr = me.ReferenceFrame('N')
# NOTE : This does not allow for rotation about the pendulum shaft's axes,
# which is a more general spherical joint.
conical_fr = newtonian_fr.orientnew('A', 'Body', (q[0], q[1], 0), 'XYZ')
plate_fr = conical_fr.orientnew('B', 'Axis', (q[3], conical_fr.z))
pendulum_fr = plate_fr.orientnew('C', 'Axis', (q[7], plate_fr.y))

## Set angular velocities.
# Here we ensure that the velocities are all expressed in terms of the
# generalized speeds.
# To set the angular velocity of the conical pendulum's quaternion
# coordinates we simply compute the angular velocity and substitute in the
# simple kinematical differential equation relations defined above.
conical_ang_vel = conical_fr.ang_vel_in(newtonian_fr)
conical_ang_vel = conical_ang_vel.subs(kin_diff_eqs)
conical_fr.set_ang_vel(newtonian_fr, conical_ang_vel)
plate_fr.set_ang_vel(conical_fr, u[3] * conical_fr.z)
pendulum_fr.set_ang_vel(plate_fr, u[7] * plate_fr.y)

# Define the positions of all the important points.
conical_base_pt = me.Point('O')
plate_mass_center = conical_base_pt.locatenew('P1', (conical_length + q[2]) * conical_fr.z)
block_mass_center = plate_mass_center.locatenew('P2', q[4] * plate_fr.x + q[5] * plate_fr.z)
pendulum_base_pt = plate_mass_center.locatenew('P3', q[6] * plate_fr.x)
pendulum_bob_pt = pendulum_base_pt.locatenew('P4', pendulum_length * pendulum_fr.z)

# Set linear velocities.
conical_base_pt.set_vel(newtonian_fr, 0)

v1= plate_mass_center.pos_from(conical_base_pt).diff(me.dynamicsymbols._t, newtonian_fr).subs(kin_diff_eqs)
plate_mass_center.set_vel(newtonian_fr, v1)
v2 = v1 + u[4] * plate_fr.x + u[5] * plate_fr.z
block_mass_center.set_vel(newtonian_fr, v2)
v3 = v2 + u[6] * plate_fr.x
pendulum_base_pt.set_vel(newtonian_fr, v3)
pendulum_bob_pt.v2pt_theory(pendulum_base_pt, newtonian_fr, pendulum_fr)

# forces

loads = []

# gravitational
loads.append((plate_mass_center, plate_mass * g * newtonian_fr.z))
loads.append((block_mass_center, block_mass * g * newtonian_fr.z))
loads.append((pendulum_base_pt, pendulum_base_mass * g * newtonian_fr.z))
loads.append((pendulum_bob_pt, pendulum_bob_mass * g * newtonian_fr.z))

# spring loads

# The plate is attached to the conical pendulum shaft by a linear and
# torsional spring.
loads.append((plate_mass_center, -k2 * q[2] * conical_fr.z))
loads.append((plate_fr, -k3 * q[3] * conical_fr.z))

# The block is attached to the plate by two springs.
loads.append((block_mass_center, -k4 * q[4] * plate_fr.x - k5 * q[5] * plate_fr.z))
loads.append((plate_mass_center, k4 * q[4] * plate_fr.x + k5 * q[5] * plate_fr.z))

# The base of the simple pendulum slides in a slot on the plate and is
# attached to the plate by a spring.
loads.append((pendulum_base_pt, -k6 * q[6] * plate_fr.x))
loads.append((plate_mass_center, k6 * q[6] * plate_fr.x))

# inertia

# plate
# I = m * (h^2 + w^2) / 12
plate_central_inertia = (me.inertia(plate_fr, plate_mass / 12 * plate_length**2,
                                              plate_mass / 6 * plate_length**2,
                                              plate_mass / 12 * plate_length**2),
                         plate_mass_center)

i_block = block_mass * block_length**2 / 6
block_central_inertia = (me.inertia(plate_fr, i_block, i_block, i_block),
                         block_mass_center)

# bodies and particles

plate = me.RigidBody('plane', plate_mass_center, plate_fr, plate_mass, plate_central_inertia)
block = me.RigidBody('block', block_mass_center, plate_fr, block_mass, block_central_inertia)
pend_base = me.Particle('pend_base', pendulum_base_pt, pendulum_base_mass)
pend_bob = me.Particle('pend_bob', pendulum_bob_pt, pendulum_bob_mass)

bodies = [plate, block, pend_base, pend_bob]

kane = me.KanesMethod(newtonian_fr, q, u, kd_eqs=[qd_i - u_i for qd_i, u_i in kin_diff_eqs.items()])
fr, frstar = kane.kanes_equations(loads, bodies)

constants = {plate_mass: 1.0,
             plate_length: 5.0,
             conical_length: 4.0,
             pendulum_length: 2.0,
             pendulum_base_mass: 1.0,
             pendulum_bob_mass: 1.0,
             block_length: 0.5,
             block_mass: 1.0,
             k2: 50.0,
             k3: 50.0,
             k4: 80.0,
             k5: 80.0,
             k6: 80.0,
             g: 9.81}

initial_conditions = {x_i: 0.0 for x_i in q + u}

initial_conditions[q[0]] = np.deg2rad(1.0)
initial_conditions[q[1]] = np.deg2rad(1.0)
initial_conditions[q[4]] = -1.0
initial_conditions[q[5]] = 1.0
initial_conditions[q[6]] = 1.0
initial_conditions[q[7]] = np.deg2rad(5.0)

time = np.linspace(0.0, 10.0, num=1000)

sys = System(kane, constants=constants,
             initial_conditions=initial_conditions, times=time)
sys.generate_ode_function(generator='cython')

frs = []
frs.append(VisualizationFrame(plate_fr.orientnew('R', 'Axis', (np.pi / 2, plate_fr.x)), plate.masscenter, Plane(length=plate_length, width=plate_length)))
frs.append(VisualizationFrame(block, Cube(block_length)))
frs.append(VisualizationFrame(newtonian_fr, pend_base, Sphere(radius=0.25)))
frs.append(VisualizationFrame(newtonian_fr, pend_bob, Sphere(radius=0.25)))

scene = Scene(newtonian_fr, conical_base_pt, *frs, system=sys)

scene.display()
