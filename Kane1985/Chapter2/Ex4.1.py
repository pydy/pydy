#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 4.1 from Kane 1985"""

from sympy.physics.mechanics import dot, dynamicsymbols, MechanicsStrPrinter
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy import solve, symbols, pi, sin, cos
from sympy.simplify.simplify import trigsimp

def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)

theta = symbols('theta:3')
x = symbols('x:3')
q = symbols('q')

A = ReferenceFrame('A')
B = A.orientnew('B', 'SPACE', theta, 'xyz')

O = Point('O')
P = O.locatenew('P', x[0] * A.x + x[1] * A.y + x[2] * A.z)
p = P.pos_from(O)

# From problem, point P is on L (span(B.x)) when:
constraint_eqs = {x[0] : q*cos(theta[1])*cos(theta[2]),
                  x[1] : q*cos(theta[1])*sin(theta[2]),
                  x[2] : -q*sin(theta[1])}

# If point P is on line L then r^{P/O} will have no components in the B.y or
# B.z directions since point O is also on line L and B.x is parallel to L.
assert(trigsimp(dot(P.pos_from(O), B.y).subs(constraint_eqs)) == 0)
assert(trigsimp(dot(P.pos_from(O), B.z).subs(constraint_eqs)) == 0)

