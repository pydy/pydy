#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercises 8.4, 8.5 from Kane 1985."""

from __future__ import division
from sympy import factor, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point, dot, dynamicsymbols
from util import partial_velocities, generalized_active_forces


## --- Declare symbols ---
# Define the system with 6 generalized speeds as follows:
# u1 = dot(v_S*_A, N.x), u2 = dot(w_S_A, N.y), u3 = dot(w_S_A, N.z)
# u4 = dot(w_C1_A, N.z), u5 = dot(w_C2_A, N.z), u6 = dot(v_S*_A, N.z)
q1, q2, q3, q4, q5, q6 = dynamicsymbols('q1:7')
q1d, q2d, q3d, q4d, q5d, q6d = dynamicsymbols('q1:7', level=1)
u1, u2, u3, u4 ,u5, u6 = dynamicsymbols('u1:7')
alpha1, alpha2, alpha3, beta1, beta2, beta3 = symbols('alpha1:4 beta1:4')
gamma1, gamma2, gamma3, delta1, delta2, delta3 = symbols('gamma1:4 delta1:4')
L, R = symbols('L R')

# --- ReferenceFrames ---
A = ReferenceFrame('A')
# N has its y-axis fixed with A and its z-axis fixed with S
N = A.orientnew('N', 'Axis', [q3, A.y])
# S does not rotate about its y-axis as its z-axis remains horizontal.
S = N.orientnew('S', 'Axis', [q4, N.z])
S.set_ang_vel(A, u2*N.y + u3*N.z)

# Axis of S coincides with axes of C1, C2.
C1 = S.orientnew('C1', 'Axis', [q1, S.z])
C2 = S.orientnew('C2', 'Axis', [q2, S.z])

# --- Define Points and set their velocities ---
pS_star = Point('S*') # mass center of shaft S
pS_star.set_vel(S, 0)
pS_star.set_vel(A, u1*N.x + u6*N.z)
# C1_star, C2_star are points on C1, C2 that are fixed in both S and C1, C2
pC1_star = pS_star.locatenew('C1*', -L*N.z)
pC1_star.set_vel(S, 0)
pC1_star.set_vel(C1, 0)
pC1_star.v2pt_theory(pS_star, A, S)
pC2_star = pS_star.locatenew('C2*', +L*N.z)
pC2_star.set_vel(S, 0)
pC2_star.set_vel(C2, 0)
pC2_star.v2pt_theory(pS_star, A, S)
# S1, S2 are points on S that are fixed in both S and C1, C2
pS1 = pC1_star.locatenew('S1^', 0)
pS1.set_vel(A, pC1_star.vel(A))
pS2 = pC2_star.locatenew('S2^', 0)
pS2.set_vel(A, pC2_star.vel(A))
# C1_hat, C2_hat are points on C1, C2 that contact plane P
pC1_hat = pC1_star.locatenew('C1^', -R*N.y)
pC1_hat.set_vel(C1, 0)
pC1_hat.v2pt_theory(pC1_star, S, C1)
pC1_hat.v2pt_theory(pC1_star, A, C1)
pC2_hat = pC2_star.locatenew('C2^', -R*N.y)
pC2_hat.set_vel(C2, 0)
pC2_hat.v2pt_theory(pC2_star, S, C2)
pC2_hat.v2pt_theory(pC2_star, A, C2)

## --- Define kinematic differential equations ---
kde = [u4 - dot(C1.ang_vel_in(A), N.z), u5 - dot(C2.ang_vel_in(A), N.z)]
kde_map = solve(kde, [q1d, q2d])

# --- velocity constraints ---
# Assume C1, C2 roll without slip so velocity at C1^ and C2^ is zero.
vc = [dot(p.vel(A).subs(kde_map), b) for b in S for p in [pC1_hat, pC2_hat]]
vc_map = solve(vc, [u4, u5, u6])

## --- Define contact forces between S and C1, C2 ---
M1 = alpha1*N.x + alpha2*N.y + alpha3*N.z # torques
M2 = beta1*N.x + beta2*N.y + beta3*N.z
K1 = gamma1*N.x + gamma2*N.y + gamma3*N.z # forces
K2 = delta1*N.x + delta2*N.y + delta3*N.z

forces = [(pS1, -K1), (pS2, -K2), (pC1_star, K1), (pC2_star, K2)]
torques = [(S, -M1 - M2), (C1, M1), (C2, M2)]
system = zip(*forces + torques)[0]

partials = partial_velocities(system, [u1, u2, u3], A, kde_map, vc_map, N)
Fr, _ = generalized_active_forces(partials, forces + torques)
print("Exercise 8.4")
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0}_tilde = {1}".format(i, factor(simplify(f))))

#For Exercise 8.5, define 2 new symbols and add 2 new velocity constraints.
omega1, omega2 = symbols('omega1 omega2')
# w_C1_S, w_C2_S are defined in terms of q1d, q2d. We want to substitute in the
# values of ui using the kinematic DEs so that u1 determined in terms of
# u2, u3, u4, u5, u6.
vc += [dot(C1.ang_vel_in(S).subs(kde_map), N.z) - omega1,
       dot(C2.ang_vel_in(S).subs(kde_map), N.z) - omega2]
vc_map = solve(vc, [u1, u2, u4, u5, u6])
partials = partial_velocities(system, [u3], A, kde_map, vc_map, N)
Fr, _ = generalized_active_forces(partials, forces + torques)

print("\nExercise 8.5")
print("Generalized active forces:")
for i, f in enumerate(Fr, 3):
    print("F{0}_tilde = {1}".format(i, factor(simplify(f))))
