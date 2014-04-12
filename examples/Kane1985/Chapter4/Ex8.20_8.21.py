#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercises 8.20, 8.21 from Kane 1985."""

from __future__ import division
from sympy import simplify, solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia, RigidBody
from sympy.physics.mechanics import cross, dot, dynamicsymbols
from util import msprint, partial_velocities, generalized_inertia_forces


# Define generalized coordinates, speeds, and constants:
q0, q1, q2 = dynamicsymbols('q0:3')
q0d, q1d, q2d = dynamicsymbols('q0:3', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')
LA, LB, LP = symbols('LA LB LP')
p1, p2, p3 = symbols('p1:4')
A1, A2, A3 = symbols('A1:4')
B1, B2, B3 = symbols('B1:4')
C1, C2, C3 = symbols('C1:4')
D11, D22, D33, D12, D23, D31 = symbols('D11 D22 D33 D12 D23 D31')
g, mA, mB, mC, mD, t = symbols('g mA mB mC mD t')
TA_star, TB_star, TC_star, TD_star = symbols('TA* TB* TC* TD*')

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

# --- define central inertias and rigid bodies ---
IA = inertia(A, A1, A2, A3)
IB = inertia(B, B1, B2, B3)
IC = inertia(B, C1, C2, C3)
ID = inertia(B, D11, D22, D33, D12, D23, D31)

# inertia[0] is defined to be the central inertia for each rigid body
rbA = RigidBody('rbA', pA_star, A, mA, (IA, pA_star))
rbB = RigidBody('rbB', pB_star, B, mB, (IB, pB_star))
rbC = RigidBody('rbC', pC_star, C, mC, (IC, pC_star))
rbD = RigidBody('rbD', pD_star, D, mD, (ID, pD_star))
bodies = [rbA, rbB, rbC, rbD]

## --- generalized speeds ---
kde = [u1 - dot(A.ang_vel_in(E), A.x),
       u2 - dot(B.ang_vel_in(A), B.y),
       u3 - dot(pC_star.vel(B), B.z)]
kde_map = solve(kde, [q0d, q1d, q2d])
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

print('\nEx8.20')
# inertia torque for a rigid body:
# T* = -dot(alpha, I) - dot(cross(omega, I), omega)
T_star = lambda rb, F: (-dot(rb.frame.ang_acc_in(F), rb.inertia[0]) -
                        dot(cross(rb.frame.ang_vel_in(F), rb.inertia[0]),
                            rb.frame.ang_vel_in(F)))
for rb in bodies:
    print('\nT* ({0}) = {1}'.format(rb, msprint(T_star(rb, E).subs(kde_map))))

print('\nEx8.21')
system = [getattr(b, i) for b in bodies for i in ['frame', 'masscenter']]
partials = partial_velocities(system, [u1, u2, u3], E, kde_map)
Fr_star, _ = generalized_inertia_forces(partials, bodies, kde_map)
for i, f in enumerate(Fr_star, 1):
    print("\nF*{0} = {1}".format(i, msprint(simplify(f))))
