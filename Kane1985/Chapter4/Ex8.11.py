#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.11 from Kane 1985.
"""

from __future__ import division
from collections import OrderedDict
from sympy import diff, simplify, symbols
from sympy import sin, cos, pi, integrate, Matrix
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols
from sympy.physics.mechanics import MechanicsStrPrinter


def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)

def subs(x, *args, **kwargs):
    if not hasattr(x, 'subs'):
        if hasattr(x, '__iter__'):
            return map(lambda x: subs(x, *args, **kwargs), x)
    return x.subs(*args, **kwargs)

def partial_velocities(system, generalized_speeds, frame,
                       kde_map=None, constraint_map=None, express_frame=None):
    partials = {}
    if express_frame is None:
        express_frame = frame

    for p in system:
        if isinstance(p, Point):
            v = p.vel(frame)
        elif isinstance(p, ReferenceFrame):
            v = p.ang_vel_in(frame)
        if kde_map is not None:
            v = v.subs(kde_map)
        if constraint_map is not None:
            v = v.subs(constraint_map)
        v_r_p = OrderedDict((u, v.diff(u, express_frame))
                            for u in generalized_speeds)
        partials[p] = v_r_p
    return partials

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
forces = [(pP, dtheta), (pP_prime, -dtheta)]
torques = []

partials = partial_velocities([pP, pP_prime], [u2, u4], A, express_frame=B)
Flist = [0, 0]
for i, u in enumerate([u2, u4]):
    integrand = 0
    for p, f in forces:
        if partials[p][u] != 0 and f != 0:
            integrand += dot(partials[p][u], f)
    Flist[i] = integrate(integrate(integrand.subs(cart_sph_map),
                                   (theta, 0, 2*pi)),
                         (phi, -pi/2, pi/2)).subs(r, R)


print("Generalized active forces:")
for f, i in zip(Flist, [2, 4]):
    print("F{0} = {1}".format(i, msprint(simplify(f)))) 
