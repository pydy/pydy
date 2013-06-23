#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.13 from Kane 1985.
"""

from __future__ import division
from sympy import simplify, symbols, trigsimp
from sympy import sin, cos, pi, integrate, Matrix
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols
from util import msprint, subs, partial_velocities, generalized_active_forces


## --- Declare symbols ---
theta = dynamicsymbols('theta')
q1, q2, q3, q4, q5, q6 = dynamicsymbols('q1:7')
q1d, q2d, q3d, q4d, q5d, q6d = dynamicsymbols('q1:7', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')
u_prime, E,  R, M, g = symbols('u\' E R M g')
x, y, z, r, theta = symbols('x y z r theta')
alpha, beta = symbols('alpha beta')

# --- Reference Frames ---
C = ReferenceFrame('C')
P = C.orientnew('P', 'axis', [theta, C.x])
P.set_ang_vel(C, u1*C.x)

## --- define points D, S*, Q on frame A and their velocities ---
pP_star = Point('P*')
pP_star.set_vel(P, 0)
pP_star.set_vel(C, u2*C.x + u3*C.y)

pQ = pP_star.locatenew('Q', x*C.x + y*C.y + z*C.z)
pQ.set_vel(P, 0)
pQ.v2pt_theory(pP_star, C, P)

## --- map from cartesian to cylindrical coordinates ---
coord_pairs = [(x, x), (y, r*cos(theta)), (z, r*sin(theta))]
coord_map = dict([(x, x),
                  (y, r*cos(theta)),
                  (z, r*sin(theta))])
J = Matrix([coord_map.values()]).jacobian([x, theta, r])
dJ = trigsimp(J.det())

## --- define contact/distance forces ---
# force for a point on ring R1, R2, R3
n = alpha + beta*cos(theta/2) # contact pressure
t = u_prime*n # kinetic friction
tau = -pQ.vel(C).subs(coord_map).normalize() # direction of friction
v = -P.y # direction of surface
point_force = sum(simplify(dot(n*v + t*tau, b)) * b for b in P)

# want to find gen. active forces for motions where u3 = 0
forces = [(pP_star, E*C.x + M*g*C.y),
          (pQ, subs(point_force, u3, 0),
           lambda i: integrate(i.subs(coord_map) * dJ,
                               (theta, -pi, pi)).subs(r, R))]
# 3 rings so repeat the last element twice more
forces += [forces[-1]] * 2
torques = []

## --- define partial velocities ---
partials = partial_velocities([f[0] for f in forces + torques], [u1, u2, u3], C)

## -- calculate generalized active forces ---
Fr, _ = generalized_active_forces(partials, forces + torques)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))

