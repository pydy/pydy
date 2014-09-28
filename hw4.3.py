#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 4 Exercise 3."""

from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import inertia, msprint
from sympy.physics.mechanics import Point, RigidBody
from sympy import pi, solve, symbols, simplify
from sympy import acos, sin, cos


# define euler angle symbols and reference frames
q1, q2 = dynamicsymbols('q1 q2')
q1d, q2d = dynamicsymbols('q1 q2', 1)
theta, r, R, m = symbols('θ r R m', real=True, positive=True)
theta_val = pi/3
q2d_val = (-(R + r*cos(theta))/r*q1d)

N = ReferenceFrame('N')
#B = N.orientnew('B', 'body', [q1, theta, q2], 'zxz')
F1 = N.orientnew('F1', 'axis', [q1, N.z])
F2 = F1.orientnew('F2', 'axis', [theta, F1.x])
B = F2.orientnew('B', 'axis', [q2, F2.z])

# velocity of the disk at the point of contact with the ground is not moving
# since the disk rolls without slipping.
pA = Point('pA') # ball bearing A
pB = pA.locatenew('pB', -R*F1.y) # ball bearing B

pA.set_vel(N, 0)
pA.set_vel(F1, 0)

pB.set_vel(F1, 0)
pB.set_vel(B, 0)
pB.v2pt_theory(pA, N, F1)

#pC.v2pt_theory(pB, N, B)
#print('\nvelocity of point C in N, v_C_N, at q1 = 0 = ')
#print(pC.vel(N).express(N).subs(q2d, q2d_val))

Ixx = m*r**2/4
Iyy = m*r**2/4
Izz = m*r**2/2
I_disc = inertia(B, Ixx, Iyy, Izz, 0, 0, 0)
rb_disc = RigidBody('Disc', pB, B, m, (I_disc, pB))

T = rb_disc.kinetic_energy(N).subs({theta: theta_val, q2d: q2d_val})
print('T = {}'.format(msprint(simplify(T))))

t = symbols('t')
dT = T.diff(symbols('t'))
print('dT/dt = {} = 0'.format(msprint(simplify(dT))))
