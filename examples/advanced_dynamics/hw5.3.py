#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 5 Exercise 3."""

from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import inertia, msprint
from sympy.physics.mechanics import Point, RigidBody
from sympy import pi, solve, symbols, simplify
from sympy import acos, sin, cos


# 2a
q1 = dynamicsymbols('q1')

px, py, pz = symbols('px py pz', real=True, positive=True)
sx, sy, sz = symbols('sx sy sz', real=True, sositive=True)
m, g, l0, k = symbols('m g l0 k', real=True, positive=True)
Ixx, Iyy, Izz = symbols('Ixx Iyy Izz', real=True, positive=True)

N = ReferenceFrame('N')
B = N.orientnew('B', 'axis', [q1, N.z])

pA = Point('A')
pA.set_vel(N, 0)

pP = pA.locatenew('P', l0*N.y - 2*l0*N.z)

pS = pP.locatenew('S', -px*B.x - pz*B.z)

I = inertia(B, Ixx, Iyy, Izz, 0, 0, 0)
rb = RigidBody('plane', pS, B, m, (I, pS))
H = rb.angular_momentum(pS, N)
print('H about S in frame N = {}'.format(msprint(H)))
print('dH/dt = {}'.format(msprint(H.dt(N))))

print('a_S_N = {}'.format(pS.acc(N)))

