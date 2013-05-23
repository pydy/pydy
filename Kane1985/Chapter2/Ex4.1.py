#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 4.1 from Kane 1985"""

from sympy.physics.mechanics import dot, dynamicsymbols, MechanicsStrPrinter
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy import solve, symbols, pi
from sympy.simplify.simplify import trigsimp

def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)

theta1, theta2, theta3 = symbols('theta1 theta2 theta3')
x1, x2, x3 = symbols('x1 x2 x3')

A = ReferenceFrame('A')
A_1 = A.orientnew('A_1', 'Axis', [theta1, A.x])
A_2 = A_1.orientnew('A_2', 'Axis', [theta2, A.y])
B = A_2.orientnew('B', 'Axis', [theta3, A.z])

O = Point('O')
P = O.locatenew('P', x1 * A.x + x2 * A.y + x3 * A.z)
p = P.pos_from(O)

# Point P is on L (span(B.x)) when:
print("{0} = 0".format(trigsimp(dot(p, B.x))))


