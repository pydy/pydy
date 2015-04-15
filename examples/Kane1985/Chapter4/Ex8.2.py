#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.2 from Kane 1985."""

from __future__ import division
from sympy import diff, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, partial_velocities, generalized_active_forces


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

print("velocities of points P1, D* in rf A:\nv_P1_A = {0}\nv_D*_A = {1}".format(
    pP1.vel(A), pDs.vel(A)))
print("velocities of points P1, D* in rf B:\nv_P1_B = {0}\nv_D*_B = {1}".format(
    pP1.vel(B), pDs.vel(B).express(E)))

# X*B.z, (Y*E.y + Z*E.z) are forces the panes of glass
# exert on P1, D* respectively
R1 = X*B.z + C*E.x - m1*g*B.y
R2 = Y*E.y + Z*E.z - C*E.x - m2*g*B.y
forces = [(pP1, R1), (pDs, R2)]
system = [f[0] for f in forces]

# solve for u1, u2, u3 in terms of q1d, q2d, q3d and substitute
kde = [u1 - dot(pP1.vel(A), E.x), u2 - dot(pP1.vel(A), E.y), u3 - q3d]
kde_map = solve(kde, [q1d, q2d, q3d])

partials = partial_velocities(system, [u1, u2, u3], A, kde_map)
Fr, _ = generalized_active_forces(partials, forces)

# use nonholonomic partial velocities to find the nonholonomic
# generalized active forces
vc = [dot(pDs.vel(B), E.y).subs(kde_map)]
vc_map = solve(vc, [u3])
partials_tilde = partial_velocities(system, [u1, u2], A, kde_map, vc_map)
Fr_tilde, _ = generalized_active_forces(partials_tilde, forces)

print("\nFor generalized speeds {0}".format(msprint(solve(kde, [u1, u2, u3]))))
print("v_r_Pi = {0}".format(msprint(partials)))
print("\nGeneralized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))
print("\nNonholonomic generalized active forces:")
for i, f in enumerate(Fr_tilde, 1):
    print("F{0}_tilde = {1}".format(i, msprint(simplify(f))))

print("\nverify results")
A31, A32 = map(lambda x: diff(vc_map[u3], x), [u1, u2])
print("F1_tilde = {0}".format(msprint(simplify(Fr[0] + A31*Fr[2]))))
print("F2_tilde = {0}".format(msprint(simplify(Fr[1] + A32*Fr[2]))))
