#!/usr/bin/env python

# This script generates the equations of motion for a double pendulum where the
# bob rotates about the pendulum rod. It can be shown to be chaotic when
# simulated.

# import sympy and the mechanics module
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import sympy.physics.mechanics as me
from pydy.system import System
from pydy.viz import Cylinder, Plane, VisualizationFrame, Scene

# Enable pretty printing.
me.init_vprinting()

# declare the constants #
# gravity
g = sym.symbols('g')
mA, mB, lB = sym.symbols('m_A, m_B, L_B')
# plate dimensions
w, h = sym.symbols('w, h')

# declare the coordinates and speeds and their derivatives #
# theta : angle of the rod
# phi : angle of the plate relative to the rod
# omega : angular speed of the rod
# alpha : angular speed of the plate
theta, phi, omega, alpha = me.dynamicsymbols('theta phi omega alpha')

# reference frames #
# create a Newtonian reference frame
N = me.ReferenceFrame('N')
# create a reference for the rod, A, and the plate, B
A = me.ReferenceFrame('A')
B = me.ReferenceFrame('B')

# orientations #
# the rod rotates with respect to the Newtonian reference frame about the x
# axis
A.orient(N, 'Axis', (theta, N.y))
# the plate rotates about the rod's primay axis
B.orient(A, 'Axis', (phi, A.z))

# positions #
# origin of the Newtonian reference frame
No = me.Point('No')
# create a point for the mass centers of the two bodies
Ao = me.Point('Ao')
Bo = me.Point('Bo')

# define the positions of the mass centers relative to the Newtonian origin
lA = (lB - h / 2) / 2
Ao.set_pos(No, lA * A.z)
Bo.set_pos(No, lB * A.z)

# kinematical differential equations #
kinDiffs = (omega - theta.diff(),
            alpha - phi.diff())

# angular velocities
A.set_ang_vel(N, omega * N.y)
B.set_ang_vel(A, alpha * A.z)

# linear velocities and accelerations #
No.set_vel(N, 0)  # the newtonian origin is fixed
Ao.v2pt_theory(No, N, A)
Bo.v2pt_theory(No, N, A)

# central inertia

IAxx = sym.S(1) / 12 * mA * (2 * lA)**2
IAyy = IAxx
IAzz = 0

IA = (me.inertia(A, IAxx, IAyy, IAzz), Ao)

IBxx = sym.S(1) / 12 * mB * h**2
IByy = sym.S(1) / 12 * mB * (w**2 + h**2)
IBzz = sym.S(1) / 12 * mB * w**2

IB = (me.inertia(B, IBxx, IByy, IBzz), Bo)

# rigid bodies
rod = me.RigidBody('rod', Ao, A, mA, IA)
plate = me.RigidBody('plate', Bo, B, mB, IB)

# forces #
# add the gravitional force to each body
rod_gravity = (Ao, mA * g * N.z)
plate_gravity = (Bo, mB * g * N.z)

# equations of motion with Kane's method

# make a tuple of the bodies and forces
bodies = (rod, plate)
loads = (rod_gravity, plate_gravity)

# create a Kane object with respect to the Newtonian reference frame
kane = me.KanesMethod(N, q_ind=(theta, phi), u_ind=(omega, alpha),
                      kd_eqs=kinDiffs)

# calculate Kane's equations
fr, frstar = kane.kanes_equations(loads, bodies)

sys = System(kane)

sys.constants = {lB: 0.2, # m
                 h: 0.1, # m
                 w: 0.2, # m
                 mA: 0.01, # kg
                 mB: 0.1, # kg
                 g: 9.81, # m/s**2
                 }

sys.initial_conditions = {theta: np.deg2rad(90.0),
                          phi: np.deg2rad(0.5),
                          omega: 0,
                          alpha: 0}

sys.times = np.linspace(0, 10, 500)

x = sys.integrate()


plt.plot(sys.times, x)
plt.legend([sym.latex(s, mode='inline') for s in sys.coordinates + sys.speeds])

# visualize

rod_shape = Cylinder(2 * lA, 0.005, color='red')
plate_shape = Plane(h, w, color='blue')

v1 = VisualizationFrame('rod',
                        A.orientnew('rod', 'Axis', (sym.pi / 2, A.x)),
                        Ao,
                        rod_shape)

v2 = VisualizationFrame('plate',
                        B.orientnew('plate', 'Body',
                                    (sym.pi / 2, sym.pi / 2, 0), 'XZX'),
                        Bo,
                        plate_shape)

scene = Scene(N, No, v1, v2, system=sys)

scene.display()
