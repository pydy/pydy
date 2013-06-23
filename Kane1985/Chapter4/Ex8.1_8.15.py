#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercises 8.1, 8.15 from Kane 1985."""

from __future__ import division
from sympy import sin, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import dot, dynamicsymbols
from util import msprint, partial_velocities
from util import generalized_active_forces, generalized_inertia_forces


g, L, m1, m2, omega, t = symbols('g L m1 m2 omega t')
C, f1, f2 = symbols('C f1 f2')
q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')

A = ReferenceFrame('A')
B = A.orientnew('B', 'Axis', [omega * t, A.y])
E = B.orientnew('E', 'Axis', [q3, B.z])

pO = Point('O')
pO.set_vel(A, 0)
pO.set_vel(B, 0)
pP1 = pO.locatenew('P1', q1 * B.x + q2 * B.y)
pP2 = pP1.locatenew('P2', L * E.x)
pP1.set_vel(E, 0)
pP1.set_vel(B, pP1.pos_from(pO).diff(t, B))
pP1.v1pt_theory(pO, A, B)
pP2.set_vel(E, 0)
pP2.v2pt_theory(pP1, A, E)

print("velocities of points P1, P2 in rf A:\nv_P1_A = {0}\nv_P2_A = {1}".format(
    pP1.vel(A), pP2.vel(A)))

# three sets of generalized speeds
u_s1 = [dot(pP1.vel(A), A.x), dot(pP1.vel(A), A.y), q3d]
u_s2 = [dot(pP1.vel(A), E.x), dot(pP1.vel(A), E.y), q3d]
u_s3 = [q1d, q2d, q3d]

# f1, f2 are forces the panes of glass exert on P1, P2 respectively
R1 = f1*B.z + C*E.x - m1*g*B.y
R2 = f2*B.z - C*E.x - m2*g*B.y

forces = [(pP1, R1), (pP2, R2)]
point_masses = [Particle('P1', pP1, m1), Particle('P2', pP2, m2)]
torques = []

ulist = [u1, u2, u3]
for uset in [u_s1, u_s2, u_s3]:
    print("\nFor generalized speeds:\n[u1, u2, u3] = {0}".format(msprint(uset)))
    # solve for u1, u2, u3 in terms of q1d, q2d, q3d and substitute
    kde = [u_i - u_expr for u_i, u_expr in zip(ulist, uset)]
    kde_map = solve(kde, [q1d, q2d, q3d])

    # include second derivatives in kde map
    for k, v in kde_map.items():
        kde_map[k.diff(t)] = v.diff(t)

    partials = partial_velocities([pP1, pP2], ulist, A, kde_map)
    Fr, _ = generalized_active_forces(partials, forces + torques)
    Fr_star, _ = generalized_inertia_forces(partials, point_masses, kde_map)
    print("Generalized active forces:")
    for i, f in enumerate(Fr, 1):
        print("F{0} = {1}".format(i, msprint(simplify(f))))
    print("Generalized inertia forces:")
    for i, f in enumerate(Fr_star, 1):
        sub_map = {}
        if uset == u_s1: # make the results easier to read
            if i == 1 or i == 3:
                sub_map = solve([u1 - u_s1[0]], [omega*q1*sin(omega*t)])
        print("F{0}* = {1}".format(i, msprint(simplify(f.subs(sub_map)))))

