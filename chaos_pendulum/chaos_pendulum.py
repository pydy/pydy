#!/usr/bin/env python

# This script generates the equations of motion for a double pendulum where the
# bob rotates about the pendulum rod. It can be shown to be chaotic when
# simulated.

# import sympy and the mechanics module
import sympy as sym
import sympy.physics.mechanics as me

# declare the constants #
# gravity
gravity = sym.symbols('g')
# center of mass length, mass and  moment of inertia of the slender rod
lA, mA, IAxx = sym.symbols('lA mA IAxx')
# center of mass length, mass and moment of inertia of the plate
lB, mB, IBxx, IByy, IBzz = sym.symbols('lB mB IBxx IByy IBzz')

## kinematics ##

# declare the coordinates and speeds and their derivatives #
# theta : angle of the rod
# phi : angle of the plate relative to the rod
# omega : angular speed of the rod
# alpha : angular speed of the plate
theta, phi, omega, alpha = me.dynamicsymbols('theta phi omega alpha')
thetad, phid, omegad, alphad = me.dynamicsymbols('theta phi omega alpha', 1)

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
B.set_ang_acc(N, B.ang_vel_in(N).dt(N))

# linear velocities and accelerations #
No.set_vel(N, 0) # the newtonian origin is fixed
Ao.set_vel(N, omega * lA * A.x)
Ao.a2pt_theory(No, N, A)
Bo.set_vel(N, omega * lB * A.x)
Bo.a2pt_theory(No, N, A)

# kinematical differential equations #
kinDiffs = [omega - thetad, alpha - phid]

## kinetics ##

# rigid bodies #
rod = me.RigidBody('rod', Ao, A, mA, (me.inertia(A, IAxx, IAxx, 0.0), Ao)) # create the empty rod object
plate = me.RigidBody('plate', Bo, B, mB, (me.inertia(B, IBxx, IByy, IBzz), Bo)) # create the empty plate object

# forces #
# add the gravitional force to each body
rodGravity = (Ao, N.z * gravity * mA)
plateGravity = (Bo, N.z * gravity * mB)

## equations of motion with Kane's method ##

# make a list of the bodies and forces
bodyList = [rod, plate]
forceList = [rodGravity, plateGravity]

# create a Kane object with respect to the Newtonian reference frame
kane = me.KanesMethod(N, q_ind=[theta, phi], u_ind=[omega, alpha], kd_eqs=kinDiffs)


# calculate Kane's equations
fr, frstar = kane.kanes_equations(forceList, bodyList)
zero = fr + frstar
# solve Kane's equations for the derivatives of the speeds
eom = sym.solvers.solve(zero, omegad, alphad)
# add the kinematical differential equations to get the equations of motion
eom.update(kane.kindiffdict())

# print the results
me.mprint(eom)
