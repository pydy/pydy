#!/usr/bin/env python

# This script generates the equations of motion for a double pendulum where the
# bob rotates about the pendulum rod. It can be shown to be chaotic when
# simulated.

import sympy as sp
import sympy.physics.mechanics as me

# declare the coordinates and speeds and their derivatives #
# theta : angle of the rod
# phi : angle of the plate
# omega : angular speed of the rod
# alpha : angular speed of the plate
theta, phi, omega, alpha = me.dynamicsymbols('theta phi omega alpha')
thetad, phid, omegad, alphad = me.dynamicsymbols('theta phi omega alpha', 1)

# declare the constants #
# gravity
gravity = sp.symbols('g')
# center of mass length, mass and  moment of inertia of the slender rod
lA, mA, IA11 = sp.symbols('lA mA IA11')
# center of mass length, mass and moment of inertia of the plate
lB, mB, IB11, IB22, IB33 = sp.symbols('lB mB IB11 IB22 IB33')

# reference frames #
# create a Newtonian reference frame
N = me.ReferenceFrame('N')
# create a reference for the rod, A, and the plate, B
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

# orientations #
# the rod rotates with respect to the Newtonian reference frame about the 2
# axis
A.orient(N, 'Axis', [theta, N.y])
# the plate rotates about the rod's primay axis
B.orient(A, 'Axis', [phi, A.z])

# positions #
# origin of the Newtonian reference frame
No = me.Point('No')
# create a point for the mass centers of the two bodies
Ao = me.Point('Ao')
Bo = me.Point('Bo')
# define the positions of the mass centers relative to the Newtonian origin
Ao.set_pos(No, lA * A.z)
Bo.set_pos(No, lB * A.z)

# angular velocities and accelerations #
A.set_ang_vel(N, omega * N.y)
B.set_ang_vel(A, alpha * A.z)

# take the derivative of the angular velocities to get angular accelerations
A.set_ang_acc(N, A.ang_vel_in(N).dt(N))
B.set_ang_acc(A, B.ang_vel_in(A).dt(A))

# linear velocities and accelerations #
No.set_vel(N, 0)
Ao.set_vel(N, omega * lA * A.x)
Ao.a2pt_theory(No, N, A)
Bo.set_vel(N, omega * lB * A.x)
Bo.a2pt_theory(No, N, A)

# kinematical differential equations #
kinDiffs = [omega - thetad, alpha - phid]

# define the inertia dyadic for each body #
IA = me.inertia(A, IA11, IA11, 0.0)
IB = me.inertia(B, IB11, IB22, IB33)

# rigid bodies #
rod = me.RigidBody() # create the empty rod object
rod.frame = A # the reference frame
rod.mass = mA # mass
rod.mc = Ao # mass center
rod.inertia = (IA, Ao) # inertia about the mass center

plate = me.RigidBody() # create the empty plate object
plate.frame = B # the reference frame
plate.mass = mB # mass
plate.mc = Bo # mass center
plate.inertia = (IB, Bo) # inertia about the mass center

# make a list of the bodies
bodyList = [rod, plate]

# forces #
# add the gravitional force to each body
forceList = [(Ao, -N.z * gravity * mA),
             (Bo, -N.z * gravity * mB)]

# equations of motion with Kane's method #
# create a Kane object with respect to the Newtonian reference frame
kane = me.Kane(N)
# set the coordinates
kane.coords([theta, phi])
# set the speeds
kane.speeds([omega, alpha])
# set the kinematical differential equations
kane.kindiffeq(kinDiffs)

# calculate Kane's equations
fr, frstar = kane.kanes_equations(forceList, bodyList)
zero = fr + frstar
# solve Kane's equations for the derivatives of the speeds
eom = sp.solvers.solve(zero, omegad, alphad)
# add the kinematical differential equations to get the equations of motion
eom.update(kane.kindiffdict())
