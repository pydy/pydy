#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.8 from Kane 1985."""

from __future__ import division
from sympy import diff, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import partial_velocities, generalized_active_forces


## --- Declare symbols ---
q1, q2 = dynamicsymbols('q1 q2')
q1d, q2d = dynamicsymbols('q1 q2', level=1)
u1, u2 = dynamicsymbols('u1 u2')
b, g, m, L, t = symbols('b g m L t')
E, I = symbols('E I')

# --- ReferenceFrames ---
N = ReferenceFrame('N')
B = N.orientnew('B', 'Axis', [-(q2 - q1)/(2*b), N.y]) # small angle approx.

# --- Define Points and set their velocities ---
pO = Point('O') # Point O is where B* would be with zero displacement.
pO.set_vel(N, 0)

# small angle approx.
pB_star = pO.locatenew('B*', -(q1 + q2)/2 * N.x)
pP1 = pO.locatenew('P1', -q1*N.x - b*N.z)
pP2 = pO.locatenew('P2', -q2*N.x + b*N.z)
for p in [pB_star, pP1, pP2]:
    p.set_vel(N, p.pos_from(pO).diff(t, N))

## --- Define kinematic differential equations/pseudo-generalized speeds ---
kde = [u1 - q1d, u2 - q2d]
kde_map = solve(kde, [q1d, q2d])

## --- Define contact/distance forces ---
M = lambda qi, qj: 12*E*I/(L**2) * (L/3 * (qj - qi)/(2*b) - qi/2)
V = lambda qi, qj: 12*E*I/(L**3) * (qi - L/2 * (qj - qi)/(2*b))

forces = [(pP1, V(q1, q2)*N.x),
          (pB_star, -m*g*N.x),
          (pP2, V(q2, q1)*N.x)]
# M2 torque is applied in the opposite direction
torques = [(B, (M(q1, q2) - M(q2, q1))*N.y)]

partials = partial_velocities([pP1, pP2, pB_star, B], [u1, u2], N, kde_map)
Fr, _ = generalized_active_forces(partials, forces + torques)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    alpha = 12*E*I/(L**3)
    k_q1 = simplify(diff(f, q1))/alpha
    k_q2 = simplify(diff(f, q2))/alpha
    k_f = simplify(f - alpha * (k_q1*q1 + k_q2*q2))
    print("F{0} = {1} * (({2})*q1 + ({3})*q2) + {4}".format(
        i, alpha, k_q1, k_q2, k_f))
