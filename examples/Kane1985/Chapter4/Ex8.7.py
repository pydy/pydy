#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.7 from Kane 1985."""

from __future__ import division
from sympy import cos, sin, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, partial_velocities, generalized_active_forces


## --- Declare symbols ---
# Define the system with 6 generalized speeds as follows:
q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')
g, m, L, t = symbols('g m L t')
Q, R, S = symbols('Q R S')

# --- ReferenceFrames ---
N = ReferenceFrame('N')

# --- Define Points and set their velocities ---
# Simplify the system to 7 points, where each point is the aggregations of
# rods that are parallel horizontally.
pO = Point('O')
pO.set_vel(N, 0)

pP1 = pO.locatenew('P1', L/2*(cos(q1)*N.x + sin(q1)*N.y))
pP2 = pP1.locatenew('P2', L/2*(cos(q1)*N.x + sin(q1)*N.y))
pP3 = pP2.locatenew('P3', L/2*(cos(q2)*N.x - sin(q2)*N.y))
pP4 = pP3.locatenew('P4', L/2*(cos(q2)*N.x - sin(q2)*N.y))
pP5 = pP4.locatenew('P5', L/2*(cos(q3)*N.x + sin(q3)*N.y))
pP6 = pP5.locatenew('P6', L/2*(cos(q3)*N.x + sin(q3)*N.y))

for p in [pP1, pP2, pP3, pP4, pP5, pP6]:
    p.set_vel(N, p.pos_from(pO).diff(t, N))

## --- Define kinematic differential equations/pseudo-generalized speeds ---
kde = [u1 - L*q1d, u2 - L*q2d, u3 - L*q3d]
kde_map = solve(kde, [q1d, q2d, q3d])

## --- Define contact/distance forces ---
forces = [(pP1, 6*m*g*N.x),
          (pP2, S*N.y + 5*m*g*N.x),
          (pP3, 6*m*g*N.x),
          (pP4, -Q*N.y + 5*m*g*N.x),
          (pP5, 6*m*g*N.x),
          (pP6, R*N.y + 5*m*g*N.x)]

partials = partial_velocities([pP1, pP2, pP3, pP4, pP5, pP6],
                              [u1, u2, u3], N, kde_map)
Fr, _ = generalized_active_forces(partials, forces)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))
