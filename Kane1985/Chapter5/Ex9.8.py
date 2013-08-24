#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.8 from Kane 1985."""

from __future__ import division
from sympy import simplify, solve, symbols, Matrix
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia, RigidBody
from sympy.physics.mechanics import cross, dot, dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, potential_energy
from util import kde_matrix, vc_matrix


q1, q2, q3, q4, q5 = dynamicsymbols('q1:6')
q1d, q2d, q3d, q4d, q5d = dynamicsymbols('q1:6', level=1)
u1, u2, u3, u4, u5 = dynamicsymbols('u1:6')

u_prime, R, M, g, e, f, theta = symbols('u\' R, M, g, e, f, theta')
a, b, mA, mB, IA, J, K, t = symbols('a b mA mB IA J K t')
IA22, IA23, IA33 = symbols('IA22 IA23 IA33')
Q1, Q2, Q3 = symbols('Q1, Q2 Q3')
TB, TC = symbols('TB TC')

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

forces = [(pS_star, -M*g*F.x), (pQ, Q1*A.x)] # no friction at point Q
torques = [(A, -TB*A.z), (A, -TC*A.z), (B, TB*A.z), (C, TC*A.z)]
partials = partial_velocities(zip(*forces + torques)[0], [u1, u2],
                              F, kde_map, vc_map, express_frame=A)
Fr, _ = generalized_active_forces(partials, forces + torques)

q = [q1, q2, q3, q4, q5]
u = [u1, u2]
n = len(q)
p = len(u)
m = n - p

if vc_map is not None:
    u += sorted(vc_map.keys(), cmp=lambda x, y: x.compare(y))

dV_dq = symbols('∂V/∂q1:{0}'.format(n + 1))
dV_eq = Matrix(Fr).T

W_sr, _ = kde_matrix(u, kde_map)
if vc_map is not None:
    A_kr, _ = vc_matrix(u, vc_map)
else:
    A_kr = Matrix.zeros(m, p)

for s in range(W_sr.shape[0]):
    dV_eq += dV_dq[s] * (W_sr[s, :p] + W_sr[s, p:]*A_kr[:, :p])

for elem in dV_eq:
    print('{0} = 0'.format(msprint(elem)))
