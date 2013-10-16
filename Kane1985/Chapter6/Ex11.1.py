#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 11.1 from Kane 1985."""

from __future__ import division
from sympy import cancel, expand, expand_trig, solve, symbols, trigsimp
from sympy import sin, cos
from sympy.physics.mechanics import ReferenceFrame, Point, Particle
from sympy.physics.mechanics import dot, dynamicsymbols, msprint
from util import generalized_active_forces, generalized_inertia_forces
from util import partial_velocities


g, L, m1, m2, omega, t = symbols('g L m1 m2 ω t')
C, f1, f2 = symbols('C f1 f2')
q1, q2, q3 = q = dynamicsymbols('q1:4')
q1d, q2d, q3d = qd = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = u = dynamicsymbols('u1:4')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'Axis', [omega * t, A.y])
E = B.orientnew('E', 'Axis', [q3, B.z])

# points and velocities
pO = Point('O')
pO.set_vel(A, 0)
pO.set_vel(B, 0)
pP1 = pO.locatenew('P1', q1*B.x + q2*B.y)
pP2 = pP1.locatenew('P2', L * E.x)
#pP1.set_vel(E, 0)
pP1.set_vel(B, pP1.pos_from(pO).diff(t, B))
pP1.v1pt_theory(pO, A, B)
#pP2.set_vel(E, 0)
pP2.v2pt_theory(pP1, A, E)

# kinematic differential equations
kde = [u1 - dot(pP1.vel(A), E.x), u2 - dot(pP1.vel(A), E.y), u3 - q3d]
kde_map = solve(kde, qd)
# include second derivatives in kde map
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

# f1, f2 are forces the panes of glass exert on P1, P2 respectively
R1 = f1*B.z + C*E.x - m1*g*B.y
R2 = f2*B.z - C*E.x - m2*g*B.y

forces = [(pP1, R1), (pP2, R2)]
system = [Particle('P1', pP1, m1), Particle('P2', pP2, m2)]

partials = partial_velocities([pP1, pP2], u, A, kde_map)
Fr, _ = generalized_active_forces(partials, forces)
Fr_star, _ = generalized_inertia_forces(partials, system, kde_map)

# dynamical equations
dyn_eq = [x + y for x, y in zip(Fr, Fr_star)]
u1d, u2d, u3d = ud = [x.diff(t) for x in u]
dyn_eq_map = solve(dyn_eq, ud)

for x in ud:
    print('{0} = {1}'.format(msprint(x),
                             msprint(cancel(trigsimp(dyn_eq_map[x])))))

u1d_expected = (-g*sin(q3) + omega**2*q1*cos(q3) + u2*u3 +
                (omega**2*cos(q3)**2 + u3**2)*L*m2/(m1 + m2))
u2d_expected = -g*cos(q3) - (omega**2*q1*sin(q3) + u3*u1)
u3d_expected = -omega**2*sin(q3)*cos(q3)
assert expand(cancel(expand_trig(dyn_eq_map[u1d] - u1d_expected))) == 0
assert expand(cancel(expand_trig(dyn_eq_map[u2d] - u2d_expected))) == 0
assert expand(expand_trig(dyn_eq_map[u3d] - u3d_expected)) == 0
