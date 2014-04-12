#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.1 from Kane 1985."""

from __future__ import division
from sympy import expand, symbols, sin, cos, trigsimp
from sympy.physics.mechanics import ReferenceFrame, Particle, Point
from sympy.physics.mechanics import dynamicsymbols, msprint


q1, q2, q3 = q = dynamicsymbols('q1:4')
q1d, q2d, q3d = qd = dynamicsymbols('q1:4', level=1)
m1, m2, L, omega, t = symbols('m1 m2 L ω t')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'axis', [omega * t, A.y])
R = B.orientnew('R', 'axis', [q3, B.z]) # reference frame of rod

# define points
pO = Point('O')
pO.set_vel(A, 0)
pP1 = pO.locatenew('P1', q1*B.x + q2*B.y)
pP2 = pP1.locatenew('P2', L*R.x)

pP1.set_vel(B, pP1.pos_from(pO).dt(B))
pP2.v2pt_theory(pP1, B, R)
pP1.v1pt_theory(pO, A, B)
pP2.v2pt_theory(pP1, A, R)

# define particles
paP1 = Particle('P1', pP1, m1)
paP2 = Particle('P2', pP2, m2)

# kinetic energy
K_S = lambda rf: sum(pa.kinetic_energy(rf) for pa in [paP1, paP2])
K_S_A = trigsimp(K_S(A))
K_S_B = trigsimp(K_S(B))
print('K_S_A = {0}'.format(msprint(K_S_A)))
print('K_S_B = {0}'.format(msprint(K_S_B)))

K_S_B_expected = ((m1 + m2)*(q1d**2 + q2d**2)/2 -
                  m2*L*(q1d*sin(q3) - q2d*cos(q3) - L*q3d/2)*q3d)
K_S_A_expected = K_S_B + omega**2/2*(m1*q1**2 + m2*(q1 + L*cos(q3))**2)

assert expand(K_S_A - K_S_A_expected) == 0
assert expand(K_S_B - K_S_B_expected) == 0
