#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.4 from Kane 1985."""

from __future__ import division
from sympy import sin, cos, diff, expand, solve, symbols, trigsimp
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, generalized_active_forces_V


g, m1, m2, L1, L2, L3, L4 = symbols('g m1 m2 L1 L2 L3 L4')
q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')

# reference frame, points, velocities
A = ReferenceFrame('A')

pO = Point('O')
pO.set_vel(A, 0)

pP1 = pO.locatenew('P1', L1 * (cos(q1)*A.x + sin(q1)*A.y))
pP1.set_vel(A, pP1.pos_from(pO).dt(A))

pP2 = pP1.locatenew('P2', L2 * (cos(q2)*A.x + sin(q2)*A.y))
pP2.set_vel(A, pP2.pos_from(pO).dt(A))

# kinematic differential equations
kde = [u1 - q1d, u2 - q2d, u3 - q3d]
kde_map = solve(kde, [q1d, q2d, q3d])

# configuration constraints
cc = [L1*cos(q1) + L2*cos(q2) - L3*cos(q3),
      L1*sin(q1) + L2*sin(q2) - L3*sin(q3) - L4]

# Differentiate configuration constraints and treat as velocity constraints.
vc = map(lambda x: diff(x, symbols('t')), cc)
vc_map = solve(subs(vc, kde_map), [u2, u3])

forces = [(pP1, m1*g*A.x), (pP2, m2*g*A.x)]
partials = partial_velocities([pP1, pP2], [u1], A, kde_map, vc_map)
Fr, _ = generalized_active_forces(partials, forces)

assert (trigsimp(expand(Fr[0])) ==
        trigsimp(expand(-g*L1*(m1*sin(q1) +
                        m2*sin(q3)*sin(q2 - q1)/sin(q2 - q3)))))

V_candidate = -g*(m1*L1*cos(q1) + m2*L3*cos(q3))
dV_dt = diff(V_candidate, symbols('t')).subs(kde_map).subs(vc_map)
Fr_ur = trigsimp(-Fr[0] * u1)
print('Show that {0} is a potential energy of the system.'.format(
        msprint(V_candidate)))
print('dV/dt = {0}'.format(msprint(dV_dt)))
print('-F1*u1 = {0}'.format(msprint(Fr_ur)))
print('dV/dt == -sum(Fr*ur, (r, 1, p)) = -F1*u1? {0}'.format(
        expand(dV_dt) == expand(Fr_ur)))

print('\nVerify Fr = {0} using V = {1}.'.format(msprint(Fr[0]),
                                                msprint(V_candidate)))
Fr_V = generalized_active_forces_V(V_candidate, [q1, q2, q3], [u1],
                                   kde_map, vc_map)
print('Fr obtained from V = {0}'.format(msprint(Fr_V)))
print('Fr == Fr_V? {0}'.format(
        trigsimp(expand(Fr[0])) == trigsimp(expand(Fr_V[0]))))
