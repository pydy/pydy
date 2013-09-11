#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercises 11.12, 11.13 from Kane 1985."""

from __future__ import division
from sympy import pi, solve, symbols, trigsimp
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import generalized_active_forces_K
from util import lagrange_equations, subs


g, m, R = symbols('g, m R')
q1, q2, q3, q4, q5 = q = dynamicsymbols('q1:6')
q1d, q2d, q3d, q4d, q5d = qd = dynamicsymbols('q1:6', level=1)
u1, u2, u3, u4, u5 = u = dynamicsymbols('u1:6')

# referenceframes
A = ReferenceFrame('A')
B_prime = A.orientnew('B_prime', 'Axis', [q1, A.z])
B = B_prime.orientnew('B', 'Axis', [pi/2 - q2, B_prime.x])
C = B.orientnew('C', 'Axis', [q3, B.z])

# points, velocities
pO = Point('O')
pO.set_vel(A, 0)

# R is the point in plane H that comes into contact with disk C.
pR = pO.locatenew('R', q4*A.x + q5*A.y)
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

# kinematic differential equations
#kde = [dot(C.ang_vel_in(A), x) - y for x, y in zip(B, u[:3])]
#kde += [x - y for x, y in zip(qd[3:], u[3:])]
#kde_map = solve(kde, qd)
kde = [x - y for x, y in zip(u, qd)]
kde_map = solve(kde, qd)
vc = map(lambda x: dot(pC_hat.vel(A), x), [A.x, A.y])
vc_map = solve(subs(vc, kde_map), [u4, u5])

# define disc rigidbody
IC = inertia(C, m*R**2/4, m*R**2/4, m*R**2/2)
rbC = RigidBody('rbC', pC_star, C, m, (IC, pC_star))
rbC.set_potential_energy(m*g*dot(pC_star.pos_from(pR), A.z))

# potential energy
V = rbC.potential_energy
print('V = {0}'.format(msprint(V)))

# kinetic energy
K = trigsimp(rbC.kinetic_energy(A).subs(kde_map).subs(vc_map))
print('K = {0}'.format(msprint(K)))

u_indep = [u1, u2, u3]
Fr = generalized_active_forces_K(K, q, u_indep, kde_map, vc_map)
# Fr + Fr* = 0 but the dynamical equations cannot be formulated by only
# kinetic energy as Fr = -Fr* for r = 1, ..., p
print('\ngeneralized active forces, Fr')
for i, x in enumerate(Fr, 1):
    print('F{0} = {1}'.format(i, msprint(x)))

L = K - V
le = lagrange_equations(L, q, u, kde_map)
print('\nLagrange\'s equations of the second kind')
for i, x in enumerate(le, 1):
    print('eq{0}: {1} = 0'.format(i, msprint(x)))
ud = map(lambda x: x.diff(symbols('t')), u)
de_map = solve(le, ud)
for k, v in de_map.iteritems():
    print('{0} = {1}'.format(msprint(k), msprint(v)))
