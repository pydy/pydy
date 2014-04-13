#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.2 from Kane 1985."""

from __future__ import division
import sys
from sympy import cos, diff, integrate, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import partial_velocities, generalized_active_forces


# Define generalized coordinates, speeds, and constants:
# q1 is the distance from the immersed end of thin rod R to the
# surface of the fluid.
# q2 is the angle between the axis of R and vertical.
q1, q2 = dynamicsymbols('q1 q2')
# A is the cross-sectional area of rod R, rho is the fluid mass density,
# g is the acceleration due to gravity
A, rho, g = symbols('A ρ g')
zeta = symbols('ζ')
C = symbols('C')

## --- reference frames ---
N = ReferenceFrame('N')
# N.z is points upward (from the immersed end to the surface).
# Let q2 be defined as a rotation about N.x.
R = N.orientnew('R', 'Axis', [q2, N.x])

## --- define buoyancy forces ---
# Assume cross-sectional area A will not significantly affect the
# displaced volume since R is a thin rod.
V = A * q1 / cos(q2)
beta = V * rho * g * N.z

# The buoyancy force acts through the center of buoyancy.
pO = Point('pO') # define point O to be at the surface of the fluid
p1 = pO.locatenew('p1', -q1 * N.z)
p2 = p1.locatenew('p2', q1 / cos(q2) / 2 * R.z)
p1.set_vel(N, p1.pos_from(pO).dt(N))
p2.v2pt_theory(p1, N, R)

forces  = [(p2, beta)]

## --- find V ---
# since qid = ui, Fr = -dV/dqr
q = [q1, q2]
partials = partial_velocities([p2], map(diff, q), N)
Fr, _ = generalized_active_forces(partials, forces)

# check if dFr/dqs = dFs/dqr for all r, s = 1, ..., n
for r in range(len(q)):
    for s in range(r + 1, len(q)):
        if Fr[r].diff(q[s]) != Fr[s].diff(q[r]):
            print('∂F{0}/∂q{1} != ∂F{1}/∂q{0}. V does not exist.'.format(r, s))
            print('∂F{0}/∂q{1} = {2}'.format(r, s, Fr[r].diff(q[s])))
            print('∂F{1}/∂q{0} = {2}'.format(r, s, Fr[s].diff(q[r])))
            sys.exit(0)

# form V using 5.1.6
alpha = symbols('α1:{0}'.format(len(q) + 1))
q_alpha = zip(q, alpha)
V = C
for i, dV_dqr in enumerate(map(lambda x: -x, Fr)):
    V += integrate(dV_dqr.subs(dict(q_alpha[i + 1:])).subs(q[i], zeta),
                   (zeta, alpha[i], q[i]))

print('V = {0}'.format(V))
## since α1, α2, C are any functions of t, set them to zero
f = list(alpha) + [C]
V = V.subs(dict(zip(f, [0] * len(f))))
print('V = {0}'.format(V))
