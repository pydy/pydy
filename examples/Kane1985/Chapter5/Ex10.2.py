#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.2 from Kane 1985."""

from __future__ import division
from sympy import collect, expand, pi, solve, symbols, trigsimp
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import subs


g, m, Px, Py, Pz, R = symbols('g m Px Py Pz R')
q = dynamicsymbols('q1:6')
qd = dynamicsymbols('q1:6', level=1)
u1, u2, u3, u4, u5 = u = dynamicsymbols('u1:6')

## --- Define ReferenceFrames ---
A = ReferenceFrame('A')
B_prime = A.orientnew('B_prime', 'Axis', [q[0], A.z])
B = B_prime.orientnew('B', 'Axis', [pi/2 - q[1], B_prime.x])
C = B.orientnew('C', 'Axis', [q[2], B.z])

## --- Define Points and their velocities ---
pO = Point('O')
pO.set_vel(A, 0)

# R is the point in plane H that comes into contact with disk C.
pR = pO.locatenew('R', q[3]*A.x + q[4]*A.y)
pR.set_vel(A, pR.pos_from(pO).dt(A))
pR.set_vel(B, 0)

# C^ is the point in disk C that comes into contact with plane H.
pC_hat = pR.locatenew('C^', 0)
pC_hat.set_vel(C, 0)

# C* is the point at the center of disk C.
pC_star = pC_hat.locatenew('C*', R*B.y)
pC_star.set_vel(C, 0)
pC_star.set_vel(B, 0)

# calculate velocities in A
pC_star.v2pt_theory(pR, A, B)
pC_hat.v2pt_theory(pC_star, A, C)

## --- Expressions for generalized speeds u1, u2, u3, u4, u5 ---
u_expr = map(lambda x: dot(C.ang_vel_in(A), x), B)
u_expr += qd[3:]
kde = [u_i - u_ex for u_i, u_ex in zip(u, u_expr)]
kde_map = solve(kde, qd)
vc = map(lambda x: dot(pC_hat.vel(A), x), [A.x, A.y])
vc_map = solve(subs(vc, kde_map), [u4, u5])

# define disc rigidbody
I_C = inertia(C, m*R**2/4, m*R**2/4, m*R**2/2)
rbC = RigidBody('rbC', pC_star, C, m, (I_C, pC_star))

# kinetic energy
K = collect(trigsimp(rbC.kinetic_energy(A).subs(kde_map).subs(vc_map)),
            m*R**2/8)
print('K = {0}'.format(msprint(K)))

K_expected = (m*R**2/8) * (5*u1**2 + u2**2 + 6*u3**2)
assert expand(K - K_expected) == 0
