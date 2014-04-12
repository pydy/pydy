#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercises 11.7, 11.8 from Kane 1985."""

from __future__ import division
from sympy import expand, solve, symbols, trigsimp, sin, cos, Matrix
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import generalized_active_forces, generalized_inertia_forces
from util import partial_velocities, subs


# 11.7
# Define generalized coordinates, speeds, and constants
q0, q1, q2 = q = dynamicsymbols('q0:3')
q0d, q1d, q2d = qd = dynamicsymbols('q0:3', level=1)
u1, u2, u3  = u = dynamicsymbols('u1:4')
LA, LB, LP = symbols('LA LB LP')
p1, p2, p3 = symbols('p1:4')
A1, A2, A3 = symbols('A1:4')
B1, B2, B3 = symbols('B1:4')
C1, C2, C3 = symbols('C1:4')
D11, D22, D33, D12, D23, D31 = symbols('D11 D22 D33 D12 D23 D31')
g, mA, mB, mC, mD, t = symbols('g mA mB mC mD t')
TA_star, TB_star, TC_star, TD_star = symbols('TA* TB* TC* TD*')

# reference frames
E = ReferenceFrame('E')
A = E.orientnew('A', 'Axis', [q0, E.x])
B = A.orientnew('B', 'Axis', [q1, A.y])
C = B.orientnew('C', 'Axis', [0, B.x])
D = C.orientnew('D', 'Axis', [0, C.x])

# points, velocities
pO = Point('O')
pA_star = pO.locatenew('A*', LA * A.z)
pP = pO.locatenew('P', LP * A.z)
pB_star = pP.locatenew('B*', LB * B.z)
pC_star = pB_star.locatenew('C*', q2 * B.z)
pD_star = pC_star.locatenew('D*', p1*B.x + p2*B.y + p3*B.z)

pO.set_vel(E, 0) # Point O is fixed in Reference Frame E
pA_star.v2pt_theory(pO, E, A) # Point A* is fixed in Reference Frame A
pP.v2pt_theory(pO, E, A) # Point P is fixed in Reference Frame A
pB_star.v2pt_theory(pP, E, B) # Point B* is fixed in Reference Frame B
# Point C* is moving in Reference Frame B
pC_star.set_vel(B, pC_star.pos_from(pB_star).diff(t, B))
pC_star.v1pt_theory(pB_star, E, B)
pD_star.set_vel(B, pC_star.vel(B)) # Point D* is fixed rel to Point C* in B
pD_star.v1pt_theory(pB_star, E, B) # Point D* is moving in Reference Frame B
# define additional points for reaction forces
pB_hat = pC_star.locatenew('B^', 0) # Point in frame B touching Point C*
pB_hat.v2pt_theory(pP, E, B)

# centrial inertias, rigid bodies
IA = inertia(A, A1, A2, A3)
IB = inertia(B, B1, B2, B3)
IC = inertia(B, C1, C2, C3)
ID = inertia(B, D11, D22, D33, D12, D23, D31)

# inertia is defined as (central inertia, mass center) for each rigid body
rbA = RigidBody('rbA', pA_star, A, mA, (IA, pA_star))
rbB = RigidBody('rbB', pB_star, B, mB, (IB, pB_star))
rbC = RigidBody('rbC', pC_star, C, mC, (IC, pC_star))
rbD = RigidBody('rbD', pD_star, D, mD, (ID, pD_star))
bodies = [rbA, rbB, rbC, rbD]

# kinematic differential equations
kde = [u1 - dot(A.ang_vel_in(E), A.x),
       u2 - dot(B.ang_vel_in(A), B.y),
       u3 - dot(pC_star.vel(B), B.z)]
kde_map = solve(kde, qd)
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

# forces, torques
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

# partial velocities
system = [x for b in bodies for x in [b.masscenter, b.frame]]
system += [f[0] for f in forces + torques]
partials = partial_velocities(system, u, E, kde_map)

# generalized active forces
Fr, _ = generalized_active_forces(partials, forces + torques)
Fr_star, _ = generalized_inertia_forces(partials, bodies, kde_map)

# dynamical equations
dyn_eq = subs([x + y for x, y in zip(Fr, Fr_star)], kde_map)
ud = [x.diff(t) for x in u]

# rewrite in the form:
# summation(X_sr * u'_r, (r, 1, 3)) = Y_s for s = 1, 2, 3
DE = Matrix(dyn_eq)
X_rs = Matrix(map(lambda x: DE.T.diff(x), ud)).T
Y_s = -expand(DE - X_rs*Matrix(ud))

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

D13 = D31
D21 = D12
k1 = B2 - B3
k2 = B3 - B1
k3 = B1 - B2
k4 = C2 - C3
k5 = C3 - C1
k6 = C1 - C2
k7 = D33 - D22
k8 = D11 - D33
k9 = D22 - D11
k10 = B1 + k1
k11 = B3 - k3
k12 = C1 + k4
k13 = C3 - k6
k14 = D11 - k7
k15 = D31 + D13
k16 = D33 + k9

Z35 = cos(q1)*B1
Z36 = Z3*k10
Z37 = Z2*Z1
Z38 = -Z37*k2
Z39 = sin(q1)*B3
Z40 = Z4*k11
Z41 = cos(q1)*C1
Z42 = Z3*k12
Z43 = -Z37*k5
Z44 = sin(q1)*C3
Z45 = Z4*k13
Z46 = Z1*Z1
Z47 = u2**2
Z48 = Z2*Z2
Z49 = D11*cos(q1) + D13*sin(q1)
Z50 = k14*Z3 + k15*Z4 - D12*Z37 + D23*(Z47 - Z48)
Z51 = D23*sin(q1) + D21*cos(q1)
Z52 = D31*(Z48 - Z46) + k8*Z37
Z53 = D33*sin(q1) + D31*cos(q1)
Z54 = k15*Z3 + k16*Z4 + D23*Z37 + D12*(Z46 - Z47)

# expected values
T1_E_A = dot(T_EA, A.x)
T2_A_B = dot(T_AB, B.y)
K3_B_C = dot(K_BC, B.z)

X11 = -(A1 + cos(q1)*(Z35 + Z41 + Z49) + sin(q1)*(Z39 + Z44 + Z53) +
        mA*LA**2 + mB*Z6**2 + mC*Z10**2 + mD*(Z13**2 + Z15**2 + Z16**2))
X12 = -(Z51 + mD*(Z13*Z14 - Z16*p1))
X13 = -mD*Z16
X21 = X12
X22 = -(B2 + C2 + D22 + mB*LB**2 + mC*Z9**2 + mD*(Z14**2 + p1**2))
X23 = mD*p1
X31 = X13
X32 = X23
X33 = -(mC + mD)
Y1 = (cos(q1)*(Z36 + Z42 + Z50) + sin(q1)*(Z40 + Z45 + Z54) + mB*Z6*Z23 +
      mC*Z10*Z27 + mD*(Z13*Z32 + Z15*Z33 + Z16*Z34) - T1_E_A)
Y2 = (Z38 + Z43 + Z52 + mB*LB*Z22 + mC*Z9*Z26 + mD*(Z14*Z32 - p1*Z34) -
      T2_A_B + g*((mB*LB + mC*Z9 + mD*Z14)*cos(q1) - mD*p1*sin(q1)))
Y3 = mC*Z28 + mD*Z34 - K3_B_C + g*(mC + mD)*sin(q1)
X_rs_expected = Matrix([[X11, X12, X13], [X21, X22, X23], [X31, X32, X33]])
Y_s_expected = Matrix([Y1, Y2, Y3])

assert trigsimp(expand(X_rs - X_rs_expected)) == Matrix.zeros(3)
assert trigsimp(expand(Y_s - Y_s_expected)) == Matrix.zeros(3, 1)

# 11.8
# If D* lies on line B*C*, then p1 = 0, p2 = 0.
# If each central principal axis of D is parallel to one of b1, b2, b3,
# then D12, D23, D31 = 0.
D_prop_map = {p1: 0, p2: 0, D12: 0, D23: 0, D31: 0}
DE_2 = Matrix(dyn_eq).subs(D_prop_map)
X_rs_2 = Matrix(map(lambda x: DE_2.T.diff(x), ud)).T
Y_s_2 = -expand(DE_2 - X_rs_2*Matrix(ud))

# Solution will have the form ud_r = Y_r/X_rr (r = 1, 2, 3)
# if the X_rs matrix is diagonal.
X_rs_expected_2 = Matrix([[X11, 0, 0],
                          [0, X22, 0],
                          [0, 0, X33]]).subs(D_prop_map)
Y_s_expected_2 = Y_s_expected.subs(D_prop_map)

assert trigsimp(expand(X_rs_2 - X_rs_expected_2)) == Matrix.zeros(3)
assert trigsimp(expand(Y_s_2 - Y_s_expected_2)) == Matrix.zeros(3, 1)
