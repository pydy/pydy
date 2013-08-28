#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 11.4 from Kane 1985."""

from __future__ import division
from sympy import expand, solve, symbols, trigsimp, cancel
from sympy import sin, cos
from sympy.physics.mechanics import ReferenceFrame, Point, RigidBody
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import generalized_active_forces, generalized_inertia_forces
from util import partial_velocities, subs


q1, q2, q3, q4, q5 = dynamicsymbols('q1:6')
q1d, q2d, q3d, q4d, q5d = dynamicsymbols('q1:6', level=1)
u1, u2, u3, u4, u5 = dynamicsymbols('u1:6')

a, b, mA, mB, IA, K, J = symbols('a b mA mB IA K J')
IA22, IA23, IA33 = symbols('IA22 IA23 IA33')
R, M, g, e, f, theta, t = symbols('R M g e f θ t')
Q1, Q2, Q3 = symbols('Q1 Q2 Q3')

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
# include second derivatives in kde map
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)
vc = [dot(p.vel(F), A.y) for p in [pB_hat, pC_hat]] + [dot(pD.vel(F), A.z)]
vc_map = solve(subs(vc, kde_map), [u3, u4, u5])

# inertias of bodies A, B, C
# IA22, IA23, IA33 are not specified in the problem statement, but are
# necessary to define an inertia object. Although the values of
# IA22, IA23, IA33 are not known in terms of the variables given in the
# problem statement, they do not appear in the general inertia terms.
inertia_A = inertia(A, IA, IA22, IA33, 0, IA23, 0)
inertia_B = inertia(B, K, K, J)
inertia_C = inertia(C, K, K, J)

# define the rigid bodies A, B, C
rbA = RigidBody('rbA', pA_star, A, mA, (inertia_A, pA_star))
rbB = RigidBody('rbB', pB_star, B, mB, (inertia_B, pB_star))
rbC = RigidBody('rbC', pC_star, C, mB, (inertia_C, pC_star))
bodies = [rbA, rbB, rbC]

forces = [(pS_star, -M*g*F.x), (pQ, Q1*A.x + Q2*A.y + Q3*A.z)]
system = ([i.masscenter for i in bodies] + [i.frame for i in bodies] +
          list(zip(*forces)[0]))
partials = partial_velocities(system, [u1, u2], F, kde_map, vc_map)

# use nonholonomic partial velocities to find the nonholonomic
# generalized active forces
Fr, _ = generalized_active_forces(partials, forces)
Fr_star, _ = generalized_inertia_forces(partials, bodies, kde_map, vc_map)

# dynamical equations
dyn_eq = subs([x + y for x, y in zip(Fr, Fr_star)], kde_map)
u1d, u2d = ud = [x.diff(t) for x in [u1, u2]]
dyn_eq_map = solve(dyn_eq, ud)

for x in ud:
    print('{0} = {1}'.format(msprint(x),
                             msprint(trigsimp(dyn_eq_map[x]))))

u1d_expected = ((f*Q3 + M*g*e*sin(theta)*cos(q1) - mA*a*u1*u2) /
                (IA + 2*J*b**2/R**2 + 2*K + mA*a**2 + 2*mB*b**2))
u2d_expected = ((Q2 + M*g*sin(theta)*sin(q1) + mA*a*u1**2) /
                (mA + 2*mB + 2*J/R**2))
assert cancel(expand(dyn_eq_map[u1d] - u1d_expected)) == 0
assert cancel(expand(dyn_eq_map[u2d] - u2d_expected)) == 0
