#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.14 from Kane 1985."""

from __future__ import division
from sympy import sin, cos, simplify, solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import dot, dynamicsymbols, msprint
from util import generalized_inertia_forces_K, subs


q1, q2, q3 = q = dynamicsymbols('q1:4')
q1d, q2d, q3d = qd = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = u = dynamicsymbols('u1:4')
L, m1, m2, omega, t = symbols('L m1 m2 ω t')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'Axis', [omega * t, A.y])
E = B.orientnew('E', 'Axis', [q3, B.z])

# points and velocities
pO = Point('O')
pO.set_vel(A, 0)
pO.set_vel(B, 0)
pP1 = pO.locatenew('P1', q1 * B.x + q2 * B.y)
pDs = pP1.locatenew('D*', L * E.x)
pP1.set_vel(E, 0)
pP1.set_vel(B, pP1.pos_from(pO).dt(B))
pP1.v1pt_theory(pO, A, B)
pDs.set_vel(E, 0)
pDs.v2pt_theory(pP1, B, E)
pDs.v2pt_theory(pP1, A, E)

# define generalized speeds and constraints
kde = [u1 - dot(pP1.vel(A), E.x), u2 - dot(pP1.vel(A), E.y), u3 - q3d]
kde_map = solve(kde, qd)
# include second derivatives in kde map
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

vc = [dot(pDs.vel(B), E.y)]
vc_map = solve(subs(vc, kde_map), [u3])

# define system of particles
system = [Particle('P1', pP1, m1), Particle('P2', pDs, m2)]

# calculate kinetic energy, generalized inertia forces
K = sum(map(lambda x: x.kinetic_energy(A), system))
Fr_tilde_star = generalized_inertia_forces_K(K, q, [u1, u2], kde_map, vc_map)

for i, f in enumerate(Fr_tilde_star, 1):
    print("F{0}* = {1}".format(i, msprint(simplify(f))))

Fr_tilde_star_expected = [((m1 + m2)*(omega**2*q1*cos(q3) - u1.diff(t)) -
                           m1*u2**2/L + m2*L*omega**2*cos(q3)**2),
                          (-m1*(u2.diff(t) + omega**2*q1*sin(q3) -
                                u1*u2/L))]
for x, y in zip(Fr_tilde_star, Fr_tilde_star_expected):
    assert simplify(x - y) == 0
