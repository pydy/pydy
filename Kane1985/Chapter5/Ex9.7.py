#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.7 from Kane 1985."""

from __future__ import division
from sympy import cos, pi, solve, trigsimp, symbols, S
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, partial_velocities, generalized_active_forces
from util import potential_energy


q1 = dynamicsymbols('q1')
q1d = dynamicsymbols('q1', level=1)
u1 = dynamicsymbols('u1')
g, m, r, theta = symbols('g m r θ')
omega, t = symbols('ω t')

# reference frames
M = ReferenceFrame('M')
# Plane containing W rotates about a vertical line through center of W.
N = M.orientnew('N', 'Axis', [omega * t, M.z])
C = N.orientnew('C', 'Axis', [q1, N.x])
R = C # R is fixed relative to C
A = C.orientnew('A', 'Axis', [-theta, R.x])
B = C.orientnew('B', 'Axis', [theta, R.x])

# points, velocities
pO = Point('O') # Point O is at the center of the circular wire
pA = pO.locatenew('A', -r * A.z)
pB = pO.locatenew('B', -r * B.z)
pR_star = pO.locatenew('R*', 1/S(2) * (pA.pos_from(pO) + pB.pos_from(pO)))

pO.set_vel(N, 0)
pO.set_vel(C, 0)
for p in [pA, pB, pR_star]:
    p.set_vel(C, 0)
    p.v1pt_theory(pO, N, C)

# kinematic differential equations
kde = [u1 - q1d]
kde_map = solve(kde, [q1d])

# contact/distance forces
forces = [(pR_star, -m*g*N.z)]

partials = partial_velocities(zip(*forces)[0], [u1], N, kde_map)
Fr, _ = generalized_active_forces(partials, forces)

V = potential_energy(Fr, [q1], [u1], kde_map)
print('V = {0}'.format(msprint(V)))
print('Setting C = 0, α1 = π/2')
V = V.subs(dict(zip(symbols('C α1'), [0, pi/2])))
print('V = {0}'.format(msprint(V)))

assert V == -m*g*r*cos(theta)*cos(q1)
