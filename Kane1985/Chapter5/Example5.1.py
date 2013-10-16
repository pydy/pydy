#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Example 5.1 from Kane 1985."""

from __future__ import division
from sympy import Dummy, Matrix
from sympy import expand, solve, symbols, trigsimp
from sympy.physics.mechanics import ReferenceFrame, Point, dot, dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, potential_energy


g, m1, m2, k, L, omega, t = symbols('g m1 m2 k L ω t')
q1, q2, q3 = dynamicsymbols('q1:4')
qd1, qd2, qd3 = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')

## --- Define ReferenceFrames ---
A = ReferenceFrame('A')
B = A.orientnew('B', 'Axis', [omega * t, A.y])
E = B.orientnew('E', 'Axis', [q3, B.z])

## --- Define Points and their velocities ---
pO = Point('O')
pO.set_vel(A, 0)

pP1 = pO.locatenew('P1', q1*B.x + q2*B.y)
pD_star = pP1.locatenew('D*', L*E.x)

pP1.set_vel(B, pP1.pos_from(pO).dt(B))
pD_star.v2pt_theory(pP1, B, E)

pP1.v1pt_theory(pO, A, B)
pD_star.v2pt_theory(pP1, A, E)

## --- Expressions for generalized speeds u1, u2, u3 ---
kde = [u1 - dot(pP1.vel(A), E.x), u2 - dot(pP1.vel(A), E.y),
       u3 - dot(E.ang_vel_in(B), E.z)]
kde_map = solve(kde, [qd1, qd2, qd3])

## --- Velocity constraints ---
vc = [dot(pD_star.vel(B), E.y)]
vc_map = solve(subs(vc, kde_map), [u3])

## --- Define forces on each point in the system ---
K = k*E.x - k/L*dot(pP1.pos_from(pO), E.y)*E.y
gravity = lambda m: -m*g*A.y
forces = [(pP1, K), (pP1, gravity(m1)), (pD_star, gravity(m2))]

## --- Calculate generalized active forces ---
partials = partial_velocities(zip(*forces)[0], [u1, u2], A,
                              kde_map, vc_map)
Fr_tilde, _ = generalized_active_forces(partials, forces)
Fr_tilde = map(expand, map(trigsimp, Fr_tilde))

print('Finding a potential energy function V.')
V = potential_energy(Fr_tilde, [q1, q2, q3], [u1, u2], kde_map, vc_map)
if V is not None:
    print('V = {0}'.format(msprint(V)))
    print('Substituting αi = 0, C = 0...')
    zero_vars = dict(zip(symbols('C α1:4'), [0] * 4))
    print('V = {0}'.format(msprint(V.subs(zero_vars))))
