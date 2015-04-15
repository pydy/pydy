#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 7.2 from Kane 1985."""

from __future__ import division
from sympy import solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import cross, dot

c1, c2, c3 = symbols('c1 c2 c3', real=True)
alpha = 10 # [N]

N = ReferenceFrame('N')
pO = Point('O')
pS = pO.locatenew('S', 10*N.x + 5*N.z) # [m]
pR = pO.locatenew('R', 10*N.x + 12*N.y) # [m]
pQ = pO.locatenew('Q', 12*N.y + 10*N.z) # [m]
pP = pO.locatenew('P', 4*N.x + 7*N.z) # [m]

A = alpha * pQ.pos_from(pP).normalize()
B = alpha * pS.pos_from(pR).normalize()
C_ = c1*N.x + c2*N.y + c3*N.z
eqs = [dot(A + B + C_, b) for b in N]
soln = solve(eqs, [c1, c2, c3])
C = sum(soln[ci] * b
        for ci, b in zip(sorted(soln, cmp=lambda x, y: x.compare(y)), N))
print("C = {0}\nA + B + C = {1}".format(C, A + B + C))

M = cross(pP.pos_from(pO), A) + cross(pR.pos_from(pO), B)
print("|M| = {0} N-m".format(M.magnitude().n()))
