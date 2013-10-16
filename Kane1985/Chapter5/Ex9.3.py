#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.3 from Kane 1985."""

from __future__ import division
from sympy import cos, diff, expand, pi, solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, potential_energy


g, m, Px, Py, Pz, R, t = symbols('g m Px Py Pz R t')
q = dynamicsymbols('q1:6')
qd = dynamicsymbols('q1:6', level=1)
u = dynamicsymbols('u1:6')

## --- Define ReferenceFrames ---
A = ReferenceFrame('A')
B_prime = A.orientnew('B_prime', 'Axis', [q[0], A.z])
B = B_prime.orientnew('B', 'Axis', [pi/2 - q[1], B_prime.x])
C = B.orientnew('C', 'Axis', [q[2], B.z])

## --- Define Points and their velocities ---
pO = Point('O')
pO.set_vel(A, 0)

# R is the point in plane H that comes into contact with disk C.
pR = pO.locatenew('R', q[3]*A.x + q[4]*A.y)
pR.set_vel(A, pR.pos_from(pO).diff(t, A))
pR.set_vel(B, 0)

# C^ is the point in disk C that comes into contact with plane H.
pC_hat = pR.locatenew('C^', 0)
pC_hat.set_vel(C, 0)

# C* is the point at the center of disk C.
pCs = pC_hat.locatenew('C*', R*B.y)
pCs.set_vel(C, 0)
pCs.set_vel(B, 0)

# calculate velocities in A
pCs.v2pt_theory(pR, A, B)
pC_hat.v2pt_theory(pCs, A, C)

## --- Expressions for generalized speeds u1, u2, u3, u4, u5 ---
u_expr = map(lambda x: dot(C.ang_vel_in(A), x), B)
u_expr += qd[3:]
kde = [u_i - u_ex for u_i, u_ex in zip(u, u_expr)]
kde_map = solve(kde, qd)

## --- Define forces on each point in the system ---
R_C_hat = Px*A.x + Py*A.y + Pz*A.z
R_Cs = -m*g*A.z
forces = [(pC_hat, R_C_hat), (pCs, R_Cs)]

## --- Calculate generalized active forces ---
partials = partial_velocities([pC_hat, pCs], u, A, kde_map)
Fr, _ = generalized_active_forces(partials, forces)

# Impose the condition that disk C is rolling without slipping
u_indep = u[:3]
u_dep = u[3:]
vc = map(lambda x: dot(pC_hat.vel(A), x), [A.x, A.y])
vc_map = solve(subs(vc, kde_map), u_dep)

partials_tilde = partial_velocities([pC_hat, pCs], u_indep, A, kde_map, vc_map)
Fr_tilde, _ = generalized_active_forces(partials_tilde, forces)
Fr_tilde = map(expand, Fr_tilde)

# solve for ∂V/∂qs using 5.1.9
V_gamma = m * g * R * cos(q[1])
print(('\nVerify V_γ = {0} is a potential energy '.format(V_gamma) +
       'contribution of γ for C.'))
V_gamma_dot = -sum(fr * ur for fr, ur in
                   zip(*generalized_active_forces(partials_tilde,
                                                  forces[1:])))
if V_gamma_dot == V_gamma.diff(t).subs(kde_map):
    print('d/dt(V_γ) == -sum(Fr_γ * ur).')
else:
    print('d/dt(V_γ) != -sum(Fr_γ * ur).')
    print('d/dt(V_γ) = {0}'.format(msprint(V_gamma.diff(t))))
    print('-sum(Fr_γ * ur) = {0}'.format(msprint(V_gamma_dot)))

#print('\nFinding a potential energy function V while C is rolling '
#      'without slip.')
#V = potential_energy(Fr_tilde, q, u_indep, kde_map, vc_map)
#if V is not None:
#    print('V = {0}'.format(V))

print('\nFinding a potential energy function V while C is rolling with slip.')
V = potential_energy(Fr, q, u, kde_map)
if V is not None:
    print('V = {0}'.format(V))

print('\nFinding a potential energy function V while C is rolling with slip '
      'without friction.')
V = potential_energy(subs(Fr, {Px: 0, Py: 0}), q, u, kde_map)
if V is not None:
    print('Define a2, C as functions of t such that the respective '
          'contributing potential terms go to zero.')
    print('V = {0}'.format(V.subs(dict(zip(symbols('C α2'), [0, pi/2])))))
