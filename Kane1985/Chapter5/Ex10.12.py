#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.12 from Kane 1985."""

from __future__ import division
from sympy import Matrix, expand, solve, symbols, trigsimp
from sympy.physics.mechanics import Point, ReferenceFrame, RigidBody
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import inertia_coefficient_matrix, partial_velocities, subs


q1, q2, q3, q4, q5 = dynamicsymbols('q1:6')
q1d, q2d, q3d, q4d, q5d = dynamicsymbols('q1:6', level=1)
u1, u2, u3, u4, u5 = dynamicsymbols('u1:6')

R, e, f, theta = symbols('R e f theta')
a, b, mA, mB, IA = symbols('a b mA mB IA')
IA22, IA23, IA33 = symbols('IA22 IA23 IA33')

# reference frames
F = ReferenceFrame('F')
P = F.orientnew('P', 'axis', [-theta, F.y])
A = P.orientnew('A', 'axis', [q1, P.x])
# define frames for wheels
B = A.orientnew('B', 'axis', [q4, A.z])
C = A.orientnew('C', 'axis', [q5, A.z])

# define points
pO = Point('O')
pO.set_vel(F, 0)
pD = pO.locatenew('D', q2*P.y + q3*P.z)
pD.set_vel(A, 0)
pD.set_vel(F, pD.pos_from(pO).dt(F))

pS_star = pD.locatenew('S*', e*A.y)
pQ = pD.locatenew('Q', f*A.y - R*A.x)
for p in [pS_star, pQ]:
    p.set_vel(A, 0)
    p.v2pt_theory(pD, F, A)

# masscenters of bodies A, B, C
pA_star = pD.locatenew('A*', a*A.y)
pB_star = pD.locatenew('B*', -b*A.z)
pC_star = pD.locatenew('C*', +b*A.z)
for p in [pA_star, pB_star, pC_star]:
    p.set_vel(A, 0)
    p.v2pt_theory(pD, F, A)

# points of B, C touching the plane P
pB_hat = pB_star.locatenew('B^', -R*A.x)
pC_hat = pC_star.locatenew('C^', -R*A.x)
pB_hat.set_vel(B, 0)
pC_hat.set_vel(C, 0)
pB_hat.v2pt_theory(pB_star, F, B)
pC_hat.v2pt_theory(pC_star, F, C)

# kinematic differential equations and velocity constraints
kde = [u1 - dot(A.ang_vel_in(F), A.x),
       u2 - dot(pD.vel(F), A.y),
       u3 - q3d,
       u4 - q4d,
       u5 - q5d]
kde_map = solve(kde, [q1d, q2d, q3d, q4d, q5d])
vc = [dot(p.vel(F), A.y) for p in [pB_hat, pC_hat]] + [dot(pD.vel(F), A.z)]
vc_map = solve(subs(vc, kde_map), [u3, u4, u5])

# inertias of bodies A, B, C
# IA22, IA23, IA33 are not specified in the problem statement, but are
# necessary to define an inertia object. Although the values of
# IA22, IA23, IA33 are not known in terms of the variables given in the
# problem statement, they do not appear in the general inertia terms.
inertia_A = inertia(A, IA, IA22, IA33, 0, IA23, 0)
K = mB*R**2/4
J = mB*R**2/2
inertia_B = inertia(B, K, K, J)
inertia_C = inertia(C, K, K, J)

# define the rigid bodies A, B, C
rbA = RigidBody('rbA', pA_star, A, mA, (inertia_A, pA_star))
rbB = RigidBody('rbB', pB_star, B, mB, (inertia_B, pB_star))
rbC = RigidBody('rbC', pC_star, C, mB, (inertia_C, pC_star))
bodies = [rbA, rbB, rbC]

system = [i.masscenter for i in bodies] + [i.frame for i in bodies]
partials = partial_velocities(system, [u1, u2], F, kde_map, vc_map)

M = trigsimp(inertia_coefficient_matrix(bodies, partials))
print('inertia_coefficient_matrix = {0}'.format(msprint(M)))

M_expected = Matrix([[IA + mA*a**2 + mB*(R**2/2 + 3*b**2), 0],
                      [0, mA + 3*mB]])
assert expand(M - M_expected) == Matrix.zeros(2)
