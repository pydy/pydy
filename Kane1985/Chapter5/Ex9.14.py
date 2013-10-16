#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.14 from Kane 1985."""

from __future__ import division
from sympy import expand, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols
from util import msprint, partial_velocities
from util import function_from_partials, generalized_active_forces


q1, q2, q3, q4, q5, q6 = q = dynamicsymbols('q1:7')
u1, u2, u3, u4, u5, u6 = u = dynamicsymbols('u1:7')
alpha, beta = symbols('α β')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'body', [q1, q2, q3], 'xyz')

# define points
pO = Point('O')
pP = pO.locatenew('P', q1*A.x + q2*A.y + q3*A.z)
pP.set_vel(A, pP.pos_from(pO).dt(A))

# kinematic differential equations
kde_map = dict(zip(map(lambda x: x.diff(), q), u))

# forces
forces = [(pP, -beta * pP.vel(A))]
torques = [(B, -alpha * B.ang_vel_in(A))]

partials_c = partial_velocities(zip(*forces + torques )[0], u, A, kde_map)
Fr_c, _ = generalized_active_forces(partials_c, forces + torques)

dissipation_function = function_from_partials(
        map(lambda x: 0 if x == 0 else -x.subs(kde_map), Fr_c),
        u,
        zero_constants=True)
from sympy import simplify, trigsimp
dissipation_function = trigsimp(dissipation_function)
#print('ℱ = {0}'.format(msprint(dissipation_function)))

omega2 = trigsimp(dot(B.ang_vel_in(A), B.ang_vel_in(A)).subs(kde_map))
v2 = trigsimp(dot(pP.vel(A), pP.vel(A)).subs(kde_map))
sym_map = dict(zip([omega2, v2], map(lambda x: x**2, symbols('ω v'))))
#print('ω**2 = {0}'.format(msprint(omega2)))
#print('v**2 = {0}'.format(msprint(v2)))
print('ℱ = {0}'.format(msprint(dissipation_function.subs(sym_map))))

dissipation_function_expected = (alpha * omega2 + beta * v2) / 2
assert expand(dissipation_function - dissipation_function_expected) == 0
