#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.10 from Kane 1985.
Answer does not match text.
"""

from __future__ import division
from sympy import solve, symbols, sin, cos, expand, trigsimp, oo
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia
from sympy.physics.mechanics import cross, dot, dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, potential_energy


q1, q2, q3, q4 = dynamicsymbols('q1:5')
q1d, q2d, q3d, q4d = dynamicsymbols('q1:5', level=1)
u1, u2, u3, u4 = dynamicsymbols('u1:5')
m, M, G, R = symbols('m M G R')
I1, I2, I3 = symbols('I1:4')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'body', [q1, q2, q3], 'xyz')

# define points
pP = Point('P')
pP.set_vel(A, 0)

pB_star = pP.locatenew('B*', q4 * A.x)
pB_star.set_vel(B, 0)
pB_star.set_vel(A, pB_star.pos_from(pP).dt(A))

# kinematic differential equations
kde = [x - y for x, y in zip([u1, u2, u3], map(B.ang_vel_in(A).dot, B))]
kde += [u4 - q4d]
kde_map = solve(kde, [q1d, q2d, q3d, q4d])

I = inertia(B, I1, I2, I3) # central inertia dyadic of B

# forces, torques due to set of gravitational forces γ
C11, C12, C13, C21, C22, C23, C31, C32, C33 = [dot(x, y)
                                               for x in A for y in B]
f = 3/M/q4**2 * ((I1*(1 - 3*C11**2) + I2*(1 - 3*C12**2) +
                  I3*(1 - 3*C13**2))/2 * A.x +
                 (I1*C21*C11 + I2*C22*C12 + I3*C23*C13) * A.y +
                 (I1*C31*C11 + I2*C32*C12 + I3*C33*C13) * A.z)
forces = [(pB_star, -G * m * M / q4**2 * (A.x + f))]
torques = [(B, cross(3 * G * m / q4**3 * A.x, dot(I, A.x)))]

partials = partial_velocities(zip(*forces + torques)[0], [u1, u2, u3, u4],
                              A, kde_map)
Fr, _ = generalized_active_forces(partials, forces + torques)

V_gamma = potential_energy(Fr, [q1, q2, q3, q4], [u1, u2, u3, u4], kde_map)
print('V_γ = {0}'.format(msprint(V_gamma.subs(q4, R))))
print('Setting C = 0, α1, α2, α3 = 0, α4 = oo')
V_gamma = V_gamma.subs(dict(zip(symbols('C α1 α2 α3 α4'), [0]*4 + [oo] )))
print('V_γ= {0}'.format(msprint(V_gamma.subs(q4, R))))

V_gamma_expected = (-3*G*m/2/R**3 * ((I1 - I3)*sin(q2)**2 +
                                     (I1 - I2)*cos(q2)**2*sin(q3)**2) +
                    G*m*M/R + G*m/2/R**3*(2*I1 - I2 + I3))

print('V_γ - V_γ_expected = {0}'.format(
        msprint(trigsimp(expand(V_gamma.subs(q4, R)) -
                         expand(V_gamma_expected)))))
assert trigsimp(expand(V_gamma.subs(q4, R) - V_gamma_expected)) == 0
