#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.15 from Kane 1985."""

from __future__ import division
from sympy import expand, sin, cos, solve, symbols, trigsimp, together
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import dot, dynamicsymbols, msprint
from util import generalized_inertia_forces, generalized_inertia_forces_K
from util import partial_velocities


q1, q2, q3 = q = dynamicsymbols('q1:4')
q1d, q2d, q3d = qd = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = u = dynamicsymbols('u1:4')
g, m, L, t = symbols('g m L t')
Q, R, S = symbols('Q R S')

# reference frame
N = ReferenceFrame('N')

# points and velocities
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

points = [pP1, pP2, pP3, pP4, pP5, pP6]
for p in points:
    p.set_vel(N, p.pos_from(pO).diff(t, N))

# kinematic differential equations
kde = [u1 - L*q1d, u2 - L*q2d, u3 - L*q3d]
kde_map = solve(kde, qd)
# include second derivatives in kde map
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

# contact/distance forces
forces = [(pP1, 6*m*g*N.x),
          (pP2, S*N.y + 5*m*g*N.x),
          (pP3, 6*m*g*N.x),
          (pP4, -Q*N.y + 5*m*g*N.x),
          (pP5, 6*m*g*N.x),
          (pP6, R*N.y + 5*m*g*N.x)]

partials = partial_velocities(points, u, N, kde_map)
system = [Particle('P{0}'.format(i), p, x*m*g)
          for i, p, x in zip(range(1, 7), points, [6, 5] * 3)]

# part a
Fr_star_a, _ = generalized_inertia_forces(partials, system, kde_map)

# part b
K = sum(map(lambda x: x.kinetic_energy(N), system))
Fr_star_b = generalized_inertia_forces_K(K, q, u, kde_map)

# part c
G = sum(P.mass * dot(P.point.acc(N), P.point.acc(N))
        for P in system).subs(kde_map) / 2
Fr_star_c = map(lambda u_r: -G.diff(u_r.diff(t)), u)


def print_fr_star(fr_star):
    for i, f in enumerate(fr_star, 1):
        print("F{0}* = {1}".format(i, msprint(trigsimp(together(f)))))


print('part a')
print_fr_star(Fr_star_a)

print('\npart b')
print_fr_star(Fr_star_b)

print('\npart c')
print_fr_star(Fr_star_c)

for x, y in zip(Fr_star_a, Fr_star_b):
    assert expand(x - y) == 0
for x, y in zip(Fr_star_a, Fr_star_c):
    assert expand(x - y) == 0
