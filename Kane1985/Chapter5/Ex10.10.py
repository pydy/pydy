#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.10 from Kane 1985."""

from __future__ import division
from sympy import expand, solve, symbols, sin, cos, S
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import generalized_active_forces, partial_velocities
from util import potential_energy


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

# kinetic energy of robot arm E
K = sum(rb.kinetic_energy(E) for rb in bodies).subs(kde_map)
print('K = {0}'.format(msprint(K)))

# find potential energy contribution of the set of gravitational forces
forces = [(pA_star, -mA*g*E.x), (pB_star, -mB*g*E.x),
          (pC_star, -mC*g*E.x), (pD_star, -mD*g*E.x)]

## --- define partial velocities ---
partials = partial_velocities([f[0] for f in forces],
                              [u1, u2, u3], E, kde_map)

## -- calculate generalized active forces ---
Fr, _ = generalized_active_forces(partials, forces)
V = potential_energy(Fr, [q0, q1, q2], [u1, u2, u3], kde_map)
#print('V = {0}'.format(msprint(V)))
print('\nSetting C = g*mD*p1, α1, α2, α3 = 0')
V = V.subs(dict(zip(symbols('C α1 α2 α3'), [g*mD*p1, 0, 0, 0] )))
print('V = {0}'.format(msprint(V)))

Z1 = u1 * cos(q1)
Z2 = u1 * sin(q1)
Z3 = -Z2 * u2
Z4 = Z1 * u2
Z5 = -LA * u1
Z6 = -(LP + LB*cos(q1))
Z7 = u2 * LB
Z8 = Z6 * u1
Z9 = LB + q2
Z10 = Z6 - q2*cos(q1)
Z11 = u2 * Z9
Z12 = Z10 * u1
Z13 = -sin(q1) * p2
Z14 = Z9 + p3
Z15 = Z10 + sin(q1)*p1 - cos(q1)*p3
Z16 = cos(q1) * p2
Z17 = Z13*u1 + Z14*u2
Z18 = Z15 * u1
Z19 = Z16*u1 - u2*p1 + u3
Z20 = u1 * Z5
Z21 = LB * sin(q1) * u2
Z22 = -Z2 * Z8
Z23 = Z21*u1 + Z2*Z7
Z24 = Z1*Z8 - u2*Z7
Z25 = Z21 - u3*cos(q1) + q2*sin(q1)*u2
Z26 = 2*u2*u3 - Z2*Z12
Z27 = Z25*u1 + Z2*Z11 - Z1*u3
Z28 = Z1*Z12 - u2*Z11
Z29 = -Z16 * u2
Z30 = Z25 + u2*(cos(q1)*p1 + sin(q1)*p3)
Z31 = Z13 * u2
Z32 = Z29*u1 + u2*(u3 + Z19) - Z2*Z18
Z33 = Z30*u1 + Z2*Z17 - Z1*Z19
Z34 = Z31*u1 + Z1*Z18 - u2*Z17

K_expected = S(1)/2*(A1*u1**2 + (B1 + C1)*Z1**2 + (B2 + C2)*u2**2 +
                     (B3 + C3)*Z2**2 + Z1*(D11*Z1 + D12*u2 + D31*Z2) +
                     u2*(D12*Z1 + D22*u2 + D23*Z2) +
                     Z2*(D31*Z1 + D23*u2 + D33*Z2) + mA*Z5**2 +
                     mB*(Z7**2 + Z8**2) + mC*(Z11**2 + Z12**2 + u3**2) +
                     mD*(Z17**2 + Z18**2 + Z19**2))
V_expected = g*((mB*LB + mC*Z9 + mD*Z14)*sin(q1) + mD*p1*cos(q1))

assert expand(K - K_expected) == 0
assert expand(V - V_expected) == 0
