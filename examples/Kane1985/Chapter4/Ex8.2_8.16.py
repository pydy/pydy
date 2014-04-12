#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercises 8.2, 8.16 from Kane 1985."""

from __future__ import division
from sympy import diff, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import dot, dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, generalized_inertia_forces


g, L, m1, m2, omega, t = symbols('g L m1 m2 omega t')
C, X, Y, Z = symbols('C X Y Z')
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
pDs = pP1.locatenew('D*', L * E.x)
pP1.set_vel(E, 0)
pP1.set_vel(B, pP1.pos_from(pO).diff(t, B))
pP1.v1pt_theory(pO, A, B)
pDs.set_vel(E, 0)
pDs.v2pt_theory(pP1, B, E)
pDs.v2pt_theory(pP1, A, E)

# expressions for generalized speeds u1, u2, u3
u_expr = [dot(pP1.vel(A), E.x), dot(pP1.vel(A), E.y), q3d]
ulist = [u1, u2, u3]

# X*B.z, (Y*E.y + Z*E.z) are forces the panes of glass
# exert on P1, D* respectively
R1 = X*B.z + C*E.x - m1*g*B.y
R2 = Y*E.y + Z*E.z - C*E.x - m2*g*B.y
resultants = [R1, R2]
forces = [(pP1, R1), (pDs, R2)]
point_masses = [Particle('P1', pP1, m1), Particle('P2', pDs, m2)]
points = [f[0] for f in forces]

# define generalized speeds
kde = [u_i - u_ex for u_i, u_ex in zip(ulist, u_expr)]
kde_map = solve(kde, [q1d, q2d, q3d])

# include second derivatives in kde map
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

# calculate partials, generalized forces
partials = partial_velocities(points, [u1, u2, u3], A, kde_map)
Fr, _ = generalized_active_forces(partials, forces)
Fr_star, _ = generalized_inertia_forces(partials, point_masses, kde_map)

# use nonholonomic partial velocities to find the nonholonomic
# generalized active forces
vc = [dot(pDs.vel(B), E.y)]
vc_map = solve(subs(vc, kde_map), [u3])
partials_tilde = partial_velocities(points, [u1, u2], A, kde_map, vc_map)
Fr_tilde, _ = generalized_active_forces(partials_tilde, forces)
Fr_tilde_star, _ = generalized_inertia_forces(partials_tilde, point_masses,
                                              kde_map, vc_map)

print("\nFor generalized speeds\n[u1, u2, u3] = {0}".format(msprint(u_expr)))
print("\nGeneralized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))
print("\nGeneralized inertia forces:")
for i, f in enumerate(Fr_star, 1):
    print("F{0}* = {1}".format(i, msprint(simplify(f))))
print("\nNonholonomic generalized active forces:")
for i, f in enumerate(Fr_tilde, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))
print("\nNonholonomic generalized inertia forces:")
for i, f in enumerate(Fr_tilde_star, 1):
    print("F{0}* = {1}".format(i, msprint(simplify(f))))

print("\nverify results")
A31, A32 = map(lambda x: diff(vc_map[u3], x), [u1, u2])
print("F1_tilde = {0}".format(msprint(simplify(Fr[0] + A31*Fr[2]))))
print("F2_tilde = {0}".format(msprint(simplify(Fr[1] + A32*Fr[2]))))
print("F1*_tilde = {0}".format(msprint(simplify(
    (Fr_star[0] + A31*Fr_star[2]).subs(vc_map)))))
print("F2*_tilde = {0}".format(msprint(simplify(
    (Fr_star[1] + A32*Fr_star[2]).subs(vc_map)))))
