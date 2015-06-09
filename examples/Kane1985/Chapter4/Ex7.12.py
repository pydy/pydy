#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 7.12 from Kane 1985."""

from __future__ import division
from sympy import solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import cross, dot


def subdict(d, keylist):
    return dict((k, d[k]) for k in keylist)


# vectors A, B have equal magnitude 10 N
alpha = 10 # [N]
beta = symbols('beta', real=True)
b1, b2, b3 = symbols('b1 b2 b3', real=True)
p1, p2, p3 = symbols('p1 p2 p3', real=True)

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

A_prime = beta * pP.pos_from(pO).normalize()
B_prime = b1*N.x + b2*N.y + b3*N.z
pB_prime = pO.locatenew("B'", p1*N.x + p2*N.y + p3*N.z)
M_prime = cross(pB_prime.pos_from(pO), B_prime)

eqs = [dot(R - A_prime - B_prime, n) for n in N]
eqs += [dot(M - M_prime, n) for n in N]

# choose pB_prime to be perpendicular to B_prime
# then pB_prime.magnitude() gives the distance d from O
# to the line of action of B'
eqs.append(dot(pB_prime.pos_from(pO), B_prime))

soln = solve(eqs, [beta, b1, b2, b3, p1, p2, p3], dict=True)[0]

print("A' = {0}".format(A_prime.subs(subdict(soln, [beta]))))
print("|A'| = {0} N".format(soln[beta].n()))
print("B' = {0}".format(B_prime.subs(subdict(soln, [b1, b2, b3]))))
print("|B'| = {0} N".format(B_prime.subs(subdict(soln,
                                               [b1, b2, b3])).magnitude().n()))
print("pB' = {0}".format(pB_prime.pos_from(pO).subs(subdict(soln,
                                                            [p1, p2, p3]))))
print("|pB'| = {0} m".format(pB_prime.pos_from(pO).subs(
            subdict(soln, [p1, p2, p3])).magnitude().n()))
