#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 5.4 from Kane 1985."""

from __future__ import division
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import cross, dot
from sympy import integrate, simplify, symbols, integrate
from sympy import sin, cos, pi


def calc_inertia_vec(rho, p, n_a, N, iv):
    integrand = rho * cross(p, cross(n_a, p))
    return sum(simplify(integrate(dot(integrand, n), iv)) * n
               for n in N)


a, b, L, l, m, h = symbols('a b L l m h', real=True, nonnegative=True)
theta = symbols('theta', real=True)
h_theta_val = {h:b*l/L, theta:2*pi*l/L}

density = m/L
N = ReferenceFrame('N')
pO = Point('O')
pP = pO.locatenew('P', h*N.x + a*cos(theta)*N.y + a*sin(theta)*N.z)

I_1 = calc_inertia_vec(density, pP.pos_from(pO).subs(h_theta_val),
                       N.x, N, (l, 0, L))
I_2 = calc_inertia_vec(density, pP.pos_from(pO).subs(h_theta_val),
                       N.y, N, (l, 0, L))
print('I_1 = {0}'.format(I_1))
print('I_2 = {0}'.format(I_2))
