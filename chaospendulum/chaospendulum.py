#!/usr/bin/env python

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
N = me.ReferenceFrame('N', indices=['1', '2', '3'])
# create a reference for the rod, A, and the plate, B
A = me.ReferenceFrame('A', indices=['1', '2', '3'])
B = me.ReferenceFrame('B', indices=['1', '2', '3'])

# orientations #
# the rod rotates with respect to the Newtonian reference frame about the 2
# axis
A.orient(N, 'Axis', [theta, N['2']])
# the plate rotates about the rod's primay axis
B.orient(A, 'Axis', [phi, N['3']])

# positions #
# origin of the Newtonian reference frame
No = me.Point('No')
# create a point for the mass centers of the two bodies
Ao = me.Point('Ao')
Bo = me.Point('Bo')
# define the positions of the mass centers relative to the Newtonian origin
Ao.set_pos(No, lA * A['3'])
Bo.set_pos(No, lB * A['3'])

# angular velocities and accelerations #
A.set_ang_vel(N, omega * N['2'])
B.set_ang_vel(A, alpha * A['3'])

A.set_ang_acc(N, A.ang_vel_in(N).dt(N))
B.set_ang_acc(A, B.ang_vel_in(A).dt(A))

# linear velocities and accelerations #
No.set_vel(N, 0)
Ao.set_vel(N, omega * lA * A['1'])
Ao.a2pt_theory(No, N, A)
Bo.set_vel(N, omega * lB * A['1'])
Bo.a2pt_theory(No, N, A)

# kinematical differential equations #
kinDiffs = [omega - thetad, alpha - phid]

# define the inertia dyadic for each body #
IA = me.inertia(A, IA11, IA11, 0.0)
IB = me.inertia(B, IB11, IB22, IB33)

# rigid bodies #
rod = me.RigidBody()
rod.mc = Ao
rod.inertia = (IA, Ao)
rod.frame = A

plate = me.RigidBody()
plate.mc = Bo
plate.inertia = (IB, Bo)
plate.frame = B

bodyList = [rod, plate]

# forces #
# add the gravitional force to each body
forceList = [(Ao, -N['3'] * gravity * mA), (Bo, -N['3'] * gravity * mB)]

# equations of motion with Kane's method #
# create a Kane object with respect to the Newtonian reference frame
kane = me.Kane(N)
# set the coordinates
kane.coords([theta, phi])
# set the speeds
kane.speeds([omega, alpha])
# set the kinematical differential equations
kane.kindiffeq(kinDiffs)

# calculate kane's equations
kane.kanes_equations(forceList, bodyList)

# solve the equations for the double dots
massMatrix = kane.mass_matrix_full
forcing = kane.forcing_full
doubledots = massMatrix.inv() * forcing
# substitute for the derivatives of the coordinates
kinDiffDict = kane.kindiffdict()
qudots = doubledots.subs(kinDiffDict)
qudots.simplify()
