#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 7.1 from Kane 1985."""

from __future__ import division
from sympy import acos, pi, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import cross, dot

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

print("A = {0},\nB = {1}".format(A, B))
print("R = {0},\nM = {1}".format(R, M))

theta = acos(dot(R, M) / (R.magnitude() * M.magnitude()))
print("theta = {0}°".format((theta * 180 / pi).n()))
