#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 5.1 from Kane 1985
"""

from __future__ import division
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot
from sympy import solve, sqrt, symbols, integrate
from sympy import sin, cos, tan, pi, acos
from sympy import Matrix, S, simplify
import scipy.integrate

def eval_v(v, N):
    return sum(dot(v, n).evalf() * n for n in N)

def integrate_v(integrand, rf, bounds):
    """Return the integral for a Vector integrand.
    """
    # integration problems if option meijerg=True is not used
    return sum(simplify(integrate(dot(integrand, n),
                                  bounds,
                                  meijerg=True)) * n for n in rf)

m, R = symbols('m R', real=True, nonnegative=True)
theta, r, x, y, z = symbols('theta r x y z', real=True)

"""Regarding Fig. P5.1 as showing two views of a body B formed by matter
distributed uniformly (a) over a surface having no planar portions and (b)
throughout a solid, determine (by integration) the coordinates x*, y*, z* of
the mass center of B.
"""

A = ReferenceFrame('A')

# part a
print("a) mass center of surface area")
SA = (2 * pi * R) * 2
rho_a = m / SA

B = A.orientnew('B', 'Axis', [theta, A.x])
p = x * A.x  + r * B.y
y_ = dot(p, A.y)
z_ = dot(p, A.z)
J = Matrix([x, y_, z_]).jacobian([x, r, theta])
dJ = simplify(J.det())
print("dJ = {0}".format(dJ))

# calculate center of mass for the cut cylinder
# ranges from x = [1, 3], y = [-1, 1], z = [-1, 1] in Fig. P5.1
mass_cc_a = rho_a * 2*pi*R
cm_cc_a = (integrate_v(integrate_v((rho_a * p * dJ).subs(r, R),
                                           A, (theta, acos(1-x),
                                               2*pi - acos(1-x))),
                               A, (x, 0, 2)) / mass_cc_a +
           A.x)
print("cm = {0}".format(cm_cc_a))

mass_cyl_a = rho_a * 2*pi*R
cm_cyl_a = A.x/S(2)

cm_a = (mass_cyl_a*cm_cyl_a + mass_cc_a*cm_cc_a) / m
print("center of mass = {0}".format(cm_a.subs(R, 1)))

# part b
print("b) mass center of volume")
V = (pi*R**2) * 2
rho_b = m / V
mass_cc_b = rho_b * pi*R**2

# calculate center of mass for the cut cylinder
# ranges from x = [1, 3], y = [-1, 1], z = [-1, 1] in Fig. P5.1
# compute the value using scipy due to issues with sympy
R_ = 1
def pos(z, y, x, coord):
    if coord == 0:
        return x
    elif coord == 1:
        return y
    elif coord == 2:
        return z
    else:
        raise ValueError
# y bounds
def ybu(x):
    return R_*(1 - x)
def ybl(x):
    return -R_
# z bounds
def zbu(x, y):
    return (R_**2 - y**2)**0.5
def zbl(x, y):
    return -1 * zbu(x, y)

cm_cc_b = 0
for i, b in enumerate(A):
    p_i = lambda z, y, x: pos(z, y, x, i)
    cm_cc_b += scipy.integrate.tplquad(p_i, 0, 2, ybl, ybu, zbl, zbu)[0] * b
cm_cc_b *= (rho_b / mass_cc_b).subs(R, R_)
cm_cc_b += A.x

#integrand = rho_b * (x*A.x + y*A.y + z*A.z)
#cm_cc_b = (integrate_v(integrate_v(integrate_v(integrand,
#                                               A, (z,
#                                                   -sqrt(R**2 - y**2),
#                                                   sqrt(R**2 - y**2))),
#                                   A, (y, -R, R*(1 - x))),
#                       A, (x, 0, 2)) / mass_cc_b +
#           A.x)
print("cm = {0}".format(cm_cc_b))

mass_cyl_b = rho_b * pi*R**2
cm_cyl_b = A.x/S(2)

cm_b = (mass_cyl_b*cm_cyl_b + mass_cc_b*cm_cc_b) / m
print("center of mass = {0}".format(eval_v(cm_b.subs(R, 1), A)))
