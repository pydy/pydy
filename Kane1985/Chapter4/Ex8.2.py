#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.2 from Kane 1985.
"""

from __future__ import division
from sympy import diff, factor, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import partial_velocity
from sympy.physics.mechanics import MechanicsStrPrinter

def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)

def subs(x, *args, **kwargs):
    if not hasattr(x, 'subs'):
        if hasattr(x, '__iter__'):
            return map(lambda x: subs(x, *args, **kwargs), x)
    return x.subs(*args, **kwargs)

def generalized_active_forces(partial_velocities, resultants):
    assert len(partial_velocities) == len(resultants)
    Flist = []
    degree = len(partial_velocities[0])
    for r in range(degree):
        f = sum(dot(v_Pi[r], R_i)
                for v_Pi, R_i in zip(partial_velocities, resultants))
        Flist.append(factor(simplify(f)))
    return Flist

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

# expressions for generalized speeds u1, u2, u3
u_expr = [dot(pP1.vel(A), E.x), dot(pP1.vel(A), E.y), q3d]
ulist = [u1, u2, u3]

# X*B.z, (Y*E.y + Z*E.z) are forces the panes of glass
# exert on P1, D* respectively
R1 = X*B.z + C*E.x - m1*g*B.y
R2 = Y*E.y + Z*E.z - C*E.x - m2*g*B.y
resultants = [R1, R2]

# solve for u1, u2, u3 in terms of q1d, q2d, q3d and substitute
kinematic_eqs = [u_i - u_ex for u_i, u_ex in zip(ulist, u_expr)]
soln = solve(kinematic_eqs, [q1d, q2d, q3d])
vlist = subs([pP1.vel(A), pDs.vel(A)], soln)

v_r_Pi = partial_velocity(vlist, ulist, A)
F1, F2, F3 = generalized_active_forces(v_r_Pi, resultants)

# use nonholonomic partial velocities to find the nonholonomic
# generalized active forces
nh_constraint = subs([dot(pDs.vel(B), E.y)], soln)
nh_constraint_soln = solve(nh_constraint, u3)
vlist_tilde = subs(vlist, nh_constraint_soln)
v_r_Pi_tilde = partial_velocity(vlist_tilde, [u1, u2], A)
F1_tilde, F2_tilde = generalized_active_forces(v_r_Pi_tilde, resultants)

print("\nFor generalized speeds [u1, u2, u3] = {0}".format(msprint(u_expr)))
print("v_P1_A = {0}".format(vlist[0]))
print("v_D*_A = {0}".format(vlist[1]))
print("v_r_Pi = {0}".format(v_r_Pi))
print("\nGeneralized active forces:")
print("F1 = {0}".format(msprint(F1)))
print("F2 = {0}".format(msprint(F2)))
print("F3 = {0}".format(msprint(F3)))
print("\nNonholonomic generalized active forces:")
print("F1_tilde = {0}".format(msprint(F1_tilde)))
print("F2_tilde = {0}".format(msprint(F2_tilde)))

print("\nverify results")
A31, A32 = map(lambda x: diff(nh_constraint_soln[u3], x), [u1, u2])
print("F1_tilde = {0}".format(msprint(F1 + A31*F3)))
print("F2_tilde = {0}".format(msprint(F2 + A32*F3)))
