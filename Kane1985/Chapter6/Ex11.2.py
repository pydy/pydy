#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 11.2 from Kane 1985."""

from __future__ import division
from sympy import expand, solve, symbols, trigsimp
from sympy import sin, cos
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import dot, dynamicsymbols, msprint
from util import generalized_active_forces, generalized_inertia_forces
from util import partial_velocities


g, L, m1, m2, omega, t = symbols('g L m1 m2 ω t')
C, X, Y, Z = symbols('C X Y Z')
q1, q2, q3 = q = dynamicsymbols('q1:4')
q1d, q2d, q3d = qd = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = u = dynamicsymbols('u1:4')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'Axis', [omega * t, A.y])
E = B.orientnew('E', 'Axis', [q3, B.z])

# points, velocities
pO = Point('O')
pO.set_vel(A, 0)
pO.set_vel(B, 0)
pP1 = pO.locatenew('P1', q1 * B.x + q2 * B.y)
pDs = pP1.locatenew('D*', L * E.x)
pP1.set_vel(E, 0)
pP1.set_vel(B, pP1.pos_from(pO).diff(t, B))
pP1.v1pt_theory(pO, A, B)
pDs.set_vel(E, 0)
pDs.v2pt_theory(pP1, B, E)
pDs.v2pt_theory(pP1, A, E)

# X*B.z, (Y*E.y + Z*E.z) are forces the panes of glass
# exert on P1, D* respectively
R1 = X*B.z + C*E.x - m1*g*B.y
R2 = Y*E.y + Z*E.z - C*E.x - m2*g*B.y
resultants = [R1, R2]
points = [pP1, pDs]
forces = [(pP1, R1), (pDs, R2)]
system = [Particle('P1', pP1, m1), Particle('P2', pDs, m2)]

# kinematic differential equations
kde = [u1 - dot(pP1.vel(A), E.x), u2 - dot(pP1.vel(A), E.y), u3 - q3d]
kde_map = solve(kde, qd)
# include second derivatives in kde map
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

# use nonholonomic partial velocities to find the nonholonomic
# generalized active forces
vc = [dot(pDs.vel(B), E.y).subs(kde_map)]
vc_map = solve(vc, [u3])
partials = partial_velocities(points, [u1, u2], A, kde_map, vc_map)
Fr, _ = generalized_active_forces(partials, forces)
Fr_star, _ = generalized_inertia_forces(partials, system, kde_map, vc_map)

# dynamical equations
dyn_eq = [x + y for x, y in zip(Fr, Fr_star)]
u1d, u2d = ud = [x.diff(t) for x in [u1, u2]]
dyn_eq_map = solve(dyn_eq, ud)

for x in ud:
    print('{0} = {1}'.format(msprint(x),
                             msprint(trigsimp(dyn_eq_map[x]))))

u1d_expected = (-g*sin(q3) + omega**2*q1*cos(q3) +
                (m2*L*omega**2*cos(q3)**2 - m1*u2**2/L)/(m1 + m2))
u2d_expected = -g*cos(q3) - omega**2*q1*sin(q3) + u1*u2/L
assert expand(trigsimp(dyn_eq_map[u1d] - u1d_expected)) == 0
assert expand(trigsimp(dyn_eq_map[u2d] - u2d_expected)) == 0
