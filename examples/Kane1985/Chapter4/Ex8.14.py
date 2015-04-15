#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.14 from Kane 1985."""

from __future__ import division
from sympy import simplify, solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols
from util import msprint, partial_velocities, generalized_active_forces


# Define generalized coordinates, speeds, and constants:
q0, q1, q2 = dynamicsymbols('q0 q1 q2')
q0d, q1d, q2d = dynamicsymbols('q0 q1 q2', level=1)
u1, u2, u3 = dynamicsymbols('u1 u2 u3')
LA, LB, LP = symbols('LA LB LP')
p1, p2, p3 = symbols('p1 p2 p3')
g, mA, mB, mC, mD, t = symbols('g mA mB mC mD t')

## --- reference frames ---
E = ReferenceFrame('E')
A = E.orientnew('A', 'Axis', [q0, E.x])
B = A.orientnew('B', 'Axis', [q1, A.y])
C = B.orientnew('C', 'Axis', [0, B.x])
D = C.orientnew('D', 'Axis', [0, C.x])

## --- points and their velocities ---
pO = Point('O')
pA_star = pO.locatenew('A*', LA * A.z)
pP = pO.locatenew('P', LP * A.z)
pB_star = pP.locatenew('B*', LB * B.z)
pC_star = pB_star.locatenew('C*', q2 * B.z)
pD_star = pC_star.locatenew('D*', p1 * B.x + p2 * B.y + p3 * B.z)

pO.set_vel(E, 0) # Point O is fixed in Reference Frame E
pA_star.v2pt_theory(pO, E, A) # Point A* is fixed in Reference Frame A
pP.v2pt_theory(pO, E, A) # Point P is fixed in Reference Frame A
pB_star.v2pt_theory(pP, E, B) # Point B* is fixed in Reference Frame B
# Point C* is moving in Reference Frame B
pC_star.set_vel(B, pC_star.pos_from(pB_star).diff(t, B))
pC_star.v1pt_theory(pB_star, E, B)
pD_star.set_vel(B, pC_star.vel(B)) # Point D* is fixed rel to Point C* in B
pD_star.v1pt_theory(pB_star, E, B) # Point D* is moving in Reference Frame B

# --- define additional points for reaction forces ---
pB_hat = pC_star.locatenew('B^', 0) # Point in frame B touching Point C*
pB_hat.v2pt_theory(pP, E, B)

## --- generalized speeds ---
kde = [u1 - dot(A.ang_vel_in(E), A.x),
       u2 - dot(B.ang_vel_in(A), B.y),
       u3 - dot(pC_star.vel(B), B.z)]
kde_map = solve(kde, [q0d, q1d, q2d])

## --- define forces, torques ---
def define_forces(c, exert_by, exert_on, express):
    return sum(x * y
               for x, y in zip(symbols('{0}_{1}/{2}_1:4'.format(
                                    c, exert_by, exert_on)),
                               express))
T_EA = define_forces('T', E, A, A)
K_EA = define_forces('K', E, A, A)
T_AB = define_forces('T', A, B, B)
K_AB = define_forces('K', A, B, B)
T_BC = define_forces('T', B, C, B)
K_BC = define_forces('K', B, C, B)

# K_AB will be applied from A onto B and -K_AB will be applied from B onto A
# at point P so these internal forces will cancel. Note point P is fixed in
# both A and B.
forces = [(pO, K_EA), (pC_star, K_BC), (pB_hat, -K_BC),
          (pA_star, -mA*g*E.x), (pB_star, -mB*g*E.x),
          (pC_star, -mC*g*E.x), (pD_star, -mD*g*E.x)]
torques = [(A, T_EA - T_AB), (B, T_AB - T_BC), (C, T_BC)]

## --- define partial velocities ---
partials = partial_velocities([f[0] for f in forces + torques],
                              [u1, u2, u3], E, kde_map)

## -- calculate generalized active forces ---
Fr, _ = generalized_active_forces(partials, forces + torques)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))
