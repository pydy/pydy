#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 7.7 from Kane 1985."""

from __future__ import division
from sympy import acos, pi, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import cross, dot

# vectors A, B have equal magnitude alpha
alpha = symbols('alpha', real=True, positive=True)

N = ReferenceFrame('N')
pO = Point('O')
pS = pO.locatenew('S', 10*N.x + 5*N.z)
pR = pO.locatenew('R', 10*N.x + 12*N.y)
pQ = pO.locatenew('Q', 12*N.y + 10*N.z)
pP = pO.locatenew('P', 4*N.x + 7*N.z)

A = alpha * pQ.pos_from(pP).normalize()
B = alpha * pS.pos_from(pR).normalize()

R = A + B
M = cross(pP.pos_from(pO), A) + cross(pR.pos_from(pO), B)
Ms = dot(R, M) * R / dot(R, R)
ps = cross(R, M) / dot(R, R)

#print("M* = {0}".format(Ms))
print("|M*|/|R| = {0}".format(Ms.magnitude() / R.magnitude()))
print("|p*| = {0}, {1}".format(ps.magnitude(), ps.magnitude().n()))
