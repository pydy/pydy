#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.11 from Kane 1985."""

from __future__ import division
from sympy import simplify, symbols
from sympy import sin, cos, pi, integrate, Matrix
from sympy.physics.mechanics import ReferenceFrame, Point, dynamicsymbols
from util import msprint, partial_velocities, generalized_active_forces


## --- Declare symbols ---
u1, u2, u3, u4, u5, u6, u7, u8, u9 = dynamicsymbols('u1:10')
c, R = symbols('c R')
x, y, z, r, phi, theta = symbols('x y z r phi theta')

# --- Reference Frames ---
A = ReferenceFrame('A')
B = ReferenceFrame('B')
C = ReferenceFrame('C')
B.set_ang_vel(A, u1 * B.x + u2 * B.y + u3 * B.z)
C.set_ang_vel(A, u4 * B.x + u5 * B.y + u6 * B.z)
C.set_ang_vel(B, C.ang_vel_in(A) - B.ang_vel_in(A))

pC_star = Point('C*')
pC_star.set_vel(C, 0)
# since radius of cavity is very small, assume C* has zero velocity in B
pC_star.set_vel(B, 0)
pC_star.set_vel(A, u7 * B.x + u8*B.y + u9*B.z)

## --- define points P, P' ---
# point on C
pP = pC_star.locatenew('P', x * B.x + y * B.y + z * B.z)
pP.set_vel(C, 0)
pP.v2pt_theory(pC_star, B, C)
pP.v2pt_theory(pC_star, A, C)
# point on B
pP_prime = pP.locatenew("P'", 0)
pP_prime.set_vel(B, 0)
pP_prime.v1pt_theory(pC_star, A, B)

## --- Define forces ---
cart_sph_map = dict([(z, r*sin(phi)),
                      (y, r*cos(phi)*sin(theta)),
                      (x, r*cos(phi)*cos(theta))])
J = Matrix([cart_sph_map.values()]).jacobian([r, phi, theta])
dJ = simplify(J.det())

dtheta = -c * pP.vel(B) * dJ
integral = lambda i: integrate(integrate(i.subs(cart_sph_map),
                                         (theta, 0, 2*pi)),
                               (phi, -pi/2, pi/2)).subs(r, R)

forces = [(pP, dtheta, integral), (pP_prime, -dtheta, integral)]
partials = partial_velocities([pP, pP_prime], [u2, u4], A, express_frame=B)
Flist, _ = generalized_active_forces(partials, forces)

print("Generalized active forces:")
for f, i in zip(Flist, [2, 4]):
    print("F{0} = {1}".format(i, msprint(simplify(f))))
