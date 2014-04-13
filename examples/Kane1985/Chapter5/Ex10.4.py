#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.4 from Kane 1985."""

from __future__ import division
from sympy import expand, expand_trig, symbols, trigsimp
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint


m, B11, B22, B33, B12, B23, B31 = symbols('m B11 B22 B33 B12 B23 B31')
q1, q2, q3, q4, q5, q6 = dynamicsymbols('q1:7')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'body', [q4, q5, q6], 'xyz')
omega = B.ang_vel_in(A)

# points B*, O
pB_star = Point('B*')
pO = pB_star.locatenew('O', q1*B.x + q2*B.y + q3*B.z)
pO.set_vel(A, 0)
pO.set_vel(B, 0)
pB_star.v2pt_theory(pO, A, B)

# rigidbody B
I_B_O = inertia(B, B11, B22, B33, B12, B23, B31)
rbB = RigidBody('rbB', pB_star, B, m, (I_B_O, pO))

# kinetic energy
K = rbB.kinetic_energy(A)
print('K = {0}'.format(msprint(trigsimp(K))))

K_expected = dot(dot(omega, I_B_O), omega)/2
print('K_expected = 1/2*omega*I*omega = {0}'.format(
        msprint(trigsimp(K_expected))))
assert expand(expand_trig(K - K_expected)) == 0
