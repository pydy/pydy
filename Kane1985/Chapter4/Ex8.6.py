#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.6 from Kane 1985.
"""

from __future__ import division
from sympy import cos, sin, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, subs, partial_velocities, generalized_active_forces


## --- Declare symbols ---
# Define the system with 6 generalized speeds as follows:
q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')
L1, L2, L3, L4 = symbols('L1:5')
g, m1, m2, t = symbols('g m1 m2 t')

# --- ReferenceFrames ---
A = ReferenceFrame('A')

# --- Define Points and set their velocities ---
pO = Point('O')
pO.set_vel(A, 0)
pP1 = pO.locatenew('P1', L1*(cos(q1)*A.x + sin(q1)*A.y))
pP1.set_vel(A, pP1.pos_from(pO).diff(t, A))
pP2 = pP1.locatenew('P2', L2*(cos(q2)*A.x + sin(q2)*A.y))
pP2.set_vel(A, pP2.pos_from(pO).diff(t, A))

## --- configuration constraints ---
cc = [L1*cos(q1) + L2*cos(q2) - L3*cos(q3),
      L1*sin(q1) + L2*sin(q2) - L3*sin(q3) - L4]

## --- Define kinematic differential equations/pseudo-generalized speeds ---
kde = [u1 - q1d, u2 - q2d, u3 - q3d]
kde_map = solve(kde, [q1d, q2d, q3d])

# --- velocity constraints ---
vc = [c.diff(t) for c in cc]
vc_map = solve(subs(vc, kde_map), [u2, u3])

## --- Define gravitational forces ---
forces = [(pP1, m1*g*A.x), (pP2, m2*g*A.x)]

partials = partial_velocities([pP1, pP2], [u1], A, kde_map, vc_map)
Fr, _ = generalized_active_forces(partials, forces)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0}_tilde = {1}".format(i, msprint(simplify(f))))

