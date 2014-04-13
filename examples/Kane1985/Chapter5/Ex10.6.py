#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.6 from Kane 1985."""

from __future__ import division
from sympy import collect, expand, expand_trig, symbols, trigsimp, sin, cos, S
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint


q1, q2, q3 = dynamicsymbols('q1, q2 q3')
#omega1, omega2, omega3 = dynamicsymbols('ω1 ω2 ω3')
q1d, q2d = dynamicsymbols('q1, q2', level=1)
m, I11, I22, I33 = symbols('m I11 I22 I33', real=True, positive=True)

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'body', [q1, q2, q3], 'xyz')

# points B*, O
pB_star = Point('B*')
pB_star.set_vel(A, 0)

# rigidbody B
I_B_Bs = inertia(B, I11, I22, I33)
rbB = RigidBody('rbB', pB_star, B, m, (I_B_Bs, pB_star))

# kinetic energy
K = rbB.kinetic_energy(A) # velocity of point B* is zero
print('K_ω = {0}'.format(msprint(K)))

print('\nSince I11, I22, I33 are the central principal moments of inertia')
print('let I_min = I11, I_max = I33')
I_min = I11
I_max = I33
H = rbB.angular_momentum(pB_star, A)
K_min = dot(H, H) / I_max / 2
K_max = dot(H, H) / I_min / 2
print('K_ω_min = {0}'.format(msprint(K_min)))
print('K_ω_max = {0}'.format(msprint(K_max)))

print('\nI11/I33, I22/I33 =< 1, since I33 >= I11, I22, so K_ω_min <= K_ω')
print('Similarly, I22/I11, I33/I11 >= 1, '
      'since I11 <= I22, I33, so K_ω_max >= K_ω')
