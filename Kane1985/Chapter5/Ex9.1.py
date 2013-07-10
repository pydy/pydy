#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.1 from Kane 1985."""

from __future__ import division
from sympy import diff, integrate, simplify, sqrt, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols
from util import partial_velocities, generalized_active_forces


# Define generalized coordinates, speeds, and constants:
q1, q2, q3 = dynamicsymbols('q1:4')
v1, v2, v3 = symbols('v1:4')
m1, m2, C, G = symbols('m1 m2 C G')
zeta = symbols('ζ')

## --- reference frames ---
N = ReferenceFrame('N')

## --- points/particles ---
r_ = q1*N.x + q2*N.y + q3*N.z
v = v1*N.x + v2*N.y + v3*N.z
p1 = Point('p1')
p1.set_vel(N, v)
p2 = p1.locatenew('p2', r_)
p2.set_vel(N, v + r_.dt(N))

## --- define gravitational forces ---
F1 = (G * m1 * m2 / dot(r_, r_))*(p2.pos_from(p1).normalize())
F2 = -F1
forces  = [(p1, F1), (p2, F2)]

## --- find V ---
# since qid = ui, Fr = -dV/dqr
q = [q1, q2, q3]
partials = partial_velocities([p1, p2], map(diff, q), N)
Fr, _ = generalized_active_forces(partials, forces)

# check if dFr/dqs = dFs/dqr for all r, s = 1, ..., n
for r in range(len(q)):
    for s in range(len(q)):
        if Fr[r].diff(q[s]) != Fr[s].diff(q[r]):
            print('∂Fr/∂qs != ∂Fs/∂qr. V does not exist.')

# form V using 5.1.6
alpha = symbols('α1:{0}'.format(len(q) + 1))
q_alpha = zip(q, alpha)
V = C
for i, dV_dqr in enumerate(map(lambda x: -x, Fr)):
    V += integrate(dV_dqr.subs(dict(q_alpha[i + 1:])).subs(q[i], zeta),
                   (zeta, alpha[i], q[i]))

r = symbols('r')
V = V.subs(sqrt(q1**2 + q2**2 + q3**2), r)
print('V = {0}'.format(simplify(V)))
# since α1, α2, α3, C are any functions of t,
# redefine C to be C + G*m1*m2/sqrt(α1**2 + α2**2 + α3**2)
V = V.subs(C + G*m1*m2/sqrt(sum(map(lambda x: x**2, alpha))), C)
print('V = {0}'.format(V))
