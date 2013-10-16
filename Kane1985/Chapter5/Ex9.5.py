#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.5 from Kane 1985."""

from __future__ import division
from sympy import sin, cos, pi, solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, subs, partial_velocities, potential_energy
from util import generalized_active_forces, generalized_active_forces_V


q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')
g, m, L, t = symbols('g m L t')

# reference frame, points, velocities
N = ReferenceFrame('N')

pO = Point('O')
pO.set_vel(N, 0)

pP1 = pO.locatenew('P1', L/2*(cos(q1)*N.x + sin(q1)*N.y))
pP2 = pP1.locatenew('P2', L/2*(cos(q1)*N.x + sin(q1)*N.y))
pP3 = pP2.locatenew('P3', L/2*(cos(q2)*N.x - sin(q2)*N.y))
pP4 = pP3.locatenew('P4', L/2*(cos(q2)*N.x - sin(q2)*N.y))
pP5 = pP4.locatenew('P5', L/2*(cos(q3)*N.x + sin(q3)*N.y))
pP6 = pP5.locatenew('P6', L/2*(cos(q3)*N.x + sin(q3)*N.y))

for p in [pP1, pP2, pP3, pP4, pP5, pP6]:
    p.set_vel(N, p.pos_from(pO).dt(N))

# kinematic differential equations
kde = [u1 - L*q1d, u2 - L*q2d, u3 - L*q3d]
kde_map = solve(kde, [q1d, q2d, q3d])

# gravity forces
forces = [(pP1, 6*m*g*N.x),
          (pP2, 5*m*g*N.x),
          (pP3, 6*m*g*N.x),
          (pP4, 5*m*g*N.x),
          (pP5, 6*m*g*N.x),
          (pP6, 5*m*g*N.x)]

# generalized active force contribution due to gravity
partials = partial_velocities(zip(*forces)[0], [u1, u2, u3], N, kde_map)
Fr, _ = generalized_active_forces(partials, forces)

print('Potential energy contribution of gravitational forces')
V = potential_energy(Fr, [q1, q2, q3], [u1, u2, u3], kde_map)
print('V = {0}'.format(msprint(V)))
print('Setting C = 0, αi = π/2')
V = V.subs(dict(zip(symbols('C α1:4'), [0] + [pi/2]*3)))
print('V = {0}\n'.format(msprint(V)))

print('Generalized active force contributions from Vγ.')
Fr_V = generalized_active_forces_V(V, [q1, q2, q3], [u1, u2, u3], kde_map)
print('Frγ = {0}'.format(msprint(Fr_V)))
print('Fr = {0}'.format(msprint(Fr)))

