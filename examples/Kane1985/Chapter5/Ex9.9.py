#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.9 from Kane 1985."""

from __future__ import division
from sympy import solve, symbols, sin, cos, expand, trigsimp
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia
from sympy.physics.mechanics import cross, dot, dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, potential_energy
from util import generalized_active_forces_V


q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')
m, M, G, R = symbols('m M G R')
I1, I2, I3 = symbols('I1:4')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'body', [q1, q2, q3], 'xyz')

# define points
pP = Point('P')
pP.set_vel(A, 0)

pB_star = pP.locatenew('B*', R * A.x)
pB_star.set_vel(B, 0)
pB_star.set_vel(A, pB_star.pos_from(pP).dt(A))

# kinematic differential equations
kde = [x - y for x, y in zip([u1, u2, u3], map(B.ang_vel_in(A).dot, B))]
kde_map = solve(kde, [q1d, q2d, q3d])

I = inertia(B, I1, I2, I3) # central inertia dyadic of B

# forces, torques due to set of gravitational forces γ
forces = [(pB_star, -G * m * M / R**2 * A.x)]
torques = [(B, cross(3 * G * m / R**3 * A.x, dot(I, A.x)))]

partials = partial_velocities(zip(*forces + torques)[0], [u1, u2, u3],
                              A, kde_map)
Fr, _ = generalized_active_forces(partials, forces + torques)

print('part a')
V_gamma = potential_energy(Fr, [q1, q2, q3], [u1, u2, u3], kde_map)
print('V_γ = {0}'.format(msprint(V_gamma)))
print('Setting C = 0, α1, α2, α3 = 0')
V_gamma = V_gamma.subs(dict(zip(symbols('C α1 α2 α3'), [0] * 4)))
print('V_γ= {0}'.format(msprint(V_gamma)))

V_gamma_expected = (-3*G*m/2/R**3 * ((I1 - I3)*sin(q2)**2 +
                                     (I1 - I2)*cos(q2)**2*sin(q3)**2))
assert expand(V_gamma) == expand(V_gamma_expected)

print('\npart b')
kde_b = [x - y for x, y in zip([u1, u2, u3], [q1d, q2d, q3d])]
kde_map_b = solve(kde_b, [q1d, q2d, q3d])
Fr_V_gamma = generalized_active_forces_V(V_gamma, [q1, q2, q3],
                                         [u1, u2, u3], kde_map_b)
for i, fr in enumerate(Fr_V_gamma, 1):
    print('(F{0})γ = {1}'.format(i, msprint(trigsimp(fr))))

Fr_V_gamma_expected = [0,
                       3*G*m/R**3 * sin(q2)*cos(q2) * (I1*cos(q3)**2 +
                                                       I2*sin(q3)**2 - I3),
                       3*G*m/R**3 * (I1 - I2) * cos(q2)**2*sin(q3)*cos(q3)]
for x, y in zip(Fr_V_gamma, Fr_V_gamma_expected):
    assert expand(trigsimp(x - y)) == 0

print('\npart c')
q4 = dynamicsymbols('q4')
q4d = dynamicsymbols('q4', level=1)
u4 = dynamicsymbols('u4')

pB_star_c = pP.locatenew('B*', q4 * A.x)
pB_star_c.set_vel(B, 0)
pB_star_c.set_vel(A, pB_star_c.pos_from(pP).dt(A))

kde_c = kde + [u4 - q4d]
kde_map_c = solve(kde_c, [q1d, q2d, q3d, q4d])

# forces, torques due to set of gravitational forces γ
forces_c = [(pB_star_c, -G * m * M / q4**2 * A.x)]
torques_c = [(B, cross(3 * G * m / q4**3 * A.x, dot(I, A.x)))]

partials_c = partial_velocities(zip(*forces_c + torques_c)[0],
                                [u1, u2, u3, u4], A, kde_map_c)
Fr_c, _ = generalized_active_forces(partials_c, forces_c + torques_c)
V_gamma_c = potential_energy(Fr_c, [q1, q2, q3, q4],
                             [u1, u2, u3, u4], kde_map_c)
assert V_gamma_c is None
