#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.13 from Kane 1985."""

from __future__ import division
from sympy import expand, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, partial_velocities
from util import function_from_partials, generalized_active_forces


q1, q2 = q = dynamicsymbols('q1:3')
q1d, q2d = qd = dynamicsymbols('q1:3', level=1)
u1, u2 = u = dynamicsymbols('u1:3')
# L' is the natural length of the springs
alpha, beta, L1, L2, k1, k2 = symbols('α β L1 L2 k1 k2',
                                      real=True, positive=True)

# reference frames
N = ReferenceFrame('N')

# define points
pO = Point('O') # point O is fixed on the wall
pB1 = pO.locatenew('B1', (L1 + q1)*N.x) # treat block 1 as a point mass
pB2 = pB1.locatenew('B2', (L2 + q2)*N.x) # treat block 2 as a point mass
pB1.set_vel(N, pB1.pos_from(pO).dt(N))
pB2.set_vel(N, pB2.pos_from(pO).dt(N))

# kinematic differential equations
kde_map = dict(zip(map(lambda x: x.diff(), q), u))

# forces
#spring_forces = [(pB1, -k1 * q1 * N.x),
#                 (pB1, k2 * q2 * N.x),
#                 (pB2, -k2 * q2 * N.x)]
dashpot_forces = [(pB1, beta * q2d * N.x),
                 (pB2, -beta * q2d * N.x),
                 (pB2, -alpha * (q1d + q2d) * N.x)]
#forces = spring_forces + dashpot_forces

partials_c = partial_velocities(zip(*dashpot_forces)[0], u, N, kde_map)
Fr_c, _ = generalized_active_forces(partials_c, dashpot_forces)
#print('generalized active forces due to dashpot forces')
#for i, fr in enumerate(Fr_c, 1):
#    print('(F{0})c = {1} = -∂ℱ/∂u{0}'.format(i, msprint(fr)))

dissipation_function = function_from_partials(
        map(lambda x: -x.subs(kde_map), Fr_c), u, zero_constants=True)
print('ℱ = {0}'.format(msprint(dissipation_function)))

dissipation_function_expected = (alpha*u1**2 + 2*alpha*u1*u2 +
                                 (alpha + beta)*u2**2)/2
assert expand(dissipation_function - dissipation_function_expected) == 0
