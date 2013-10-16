#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.5 from Kane 1985.
Answer does not match text.
"""

from __future__ import division
from sympy import expand, symbols, trigsimp, sin, cos, S
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint


m, a, b = symbols('m a b')
q1, q2 = dynamicsymbols('q1, q2')
q1d, q2d = dynamicsymbols('q1, q2', level=1)

# reference frames
# N.x parallel to horizontal line, N.y parallel to line AC
N = ReferenceFrame('N')
A = N.orientnew('A', 'axis', [-q1, N.y])
B = A.orientnew('B', 'axis', [-q2, A.x])

# points B*, O
pO = Point('O')
pB_star = pO.locatenew('B*', S(1)/3*(2*a*B.x - b*B.y))
pO.set_vel(N, 0)
pB_star.v2pt_theory(pO, N, B)

# rigidbody B
I_B_Bs = inertia(B, m*b**2/18, m*a**2/18, m*(a**2 + b**2)/18)
rbB = RigidBody('rbB', pB_star, B, m, (I_B_Bs, pB_star))

# kinetic energy
K = rbB.kinetic_energy(N)
print('K = {0}'.format(msprint(trigsimp(K))))

K_expected = m/4*((a**2 + b**2*sin(q2)**2/3)*q1d**2 +
                  a*b*cos(q2)*q1d*q2d + b**2*q2d**2/3)
print('diff = {0}'.format(msprint(expand(trigsimp(K - K_expected)))))
assert expand(trigsimp(K - K_expected)) == 0
