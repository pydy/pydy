#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 3 Exercise 3."""

from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import msprint
from sympy.physics.mechanics import Point
from sympy import pi, solve, symbols, simplify
from sympy import acos, sin, sqrt
from sympy import Matrix

#
N = ReferenceFrame('N')

# frame B
u = (N.x + 2*N.y + -2*N.z)/3
v = (N.y + N.z)/sqrt(2)
w = (4*N.x + -1*N.y + N.z)/(3*sqrt(2))

# frame F
r = (3*N.x + -2*N.y + sqrt(3)*N.z)/4
s = (-1*N.x + sqrt(3)*N.z)/2
t = -(sqrt(3)*N.x + 2*sqrt(3)*N.y + N.z)/4

R_B_F = Matrix([[u&r, u&s, u&t],
                [v&r, v&s, v&t],
                [w&r, w&s, w&r]])
print('')
print(R_B_F)

print('')
print(R_B_F.n())
