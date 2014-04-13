#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.8 from Kane 1985."""

from __future__ import division
from sympy import sin, cos, pi, expand, simplify, solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, subs, partial_velocities, potential_energy
from util import generalized_active_forces, generalized_active_forces_V


q1, q2 = dynamicsymbols('q1 q2')
q1d, q2d = dynamicsymbols('q1 q2', level=1)
u1, u2 = dynamicsymbols('u1 u2')
b, g, m, L, t = symbols('b g m L t')
E, I = symbols('E I')

# reference frames, points, velocities
N = ReferenceFrame('N')
B = N.orientnew('B', 'Axis', [-(q2 - q1)/(2*b), N.y]) # small angle approx.

pO = Point('O') # Point O is where B* would be with zero displacement.
pO.set_vel(N, 0)

# small angle approx.
pB_star = pO.locatenew('B*', -(q1 + q2)/2 * N.x)
pP1 = pO.locatenew('P1', -q1*N.x - b*N.z)
pP2 = pO.locatenew('P2', -q2*N.x + b*N.z)
for p in [pB_star, pP1, pP2]:
    p.set_vel(N, p.pos_from(pO).diff(t, N))

# kinematic differential equations
kde = [u1 - q1d, u2 - q2d]
kde_map = solve(kde, [q1d, q2d])

# contact/distance forces
M = lambda qi, qj: 12*E*I/(L**2) * (L/3 * (qj - qi)/(2*b) - qi/2)
V = lambda qi, qj: 12*E*I/(L**3) * (qi - L/2 * (qj - qi)/(2*b))

forces = [(pP1, V(q1, q2)*N.x),
          (pB_star, -m*g*N.x),
          (pP2, V(q2, q1)*N.x)]
# M2 torque is applied in the opposite direction
torques = [(B, (M(q1, q2) - M(q2, q1))*N.y)]

partials = partial_velocities([pP1, pP2, pB_star, B], [u1, u2], N, kde_map)
Fr, _ = generalized_active_forces(partials, forces + torques)
V = simplify(potential_energy(Fr, [q1, q2], [u1, u2], kde_map))
print('V = {0}'.format(msprint(V)))
print('Setting C = 0, αi = 0')
V = V.subs(dict(zip(symbols('C α1:3'), [0] * 3)))
print('V = {0}\n'.format(msprint(V)))

assert (expand(V) ==
        expand(6*E*I/L**3 * ((1 + L/2/b + L**2/6/b**2)*(q1**2 + q2**2) -
                             q1*q2*L/b * (1 + L/3/b)) -
               m*g/2 * (q1 + q2)))
