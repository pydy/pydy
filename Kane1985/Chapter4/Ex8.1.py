#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.1 from Kane 1985."""

from __future__ import division
from sympy import factor, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols, partial_velocity
from util import msprint, subs


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

ulist = [u1, u2, u3]
for uset in [u_s1, u_s2, u_s3]:
    # solve for u1, u2, u3 in terms of q1d, q2d, q3d and substitute
    kinematic_eqs = [u_i - u_expr for u_i, u_expr in zip(ulist, uset)]
    soln = solve(kinematic_eqs, [q1d, q2d, q3d])
    vlist = subs([pP1.vel(A), pP2.vel(A)], soln)

    v_r_Pi = partial_velocity(vlist, ulist, A)
    F1, F2, F3 = [simplify(factor(
        sum(dot(v_Pi[r], R_i) for v_Pi, R_i in zip(v_r_Pi, [R1, R2]))))
        for r in range(3)]

    print("\nFor generalized speeds [u1, u2, u3] = {0}".format(msprint(uset)))
    print("v_P1_A = {0}".format(vlist[0]))
    print("v_P2_A = {0}".format(vlist[1]))
    print("v_r_Pi = {0}".format(v_r_Pi))
    print("F1 = {0}".format(msprint(F1)))
    print("F2 = {0}".format(msprint(F2)))
    print("F3 = {0}".format(msprint(F3)))
