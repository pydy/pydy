#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 7 Exercise 4."""

from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import inertia, msprint
from sympy.physics.mechanics import Point, RigidBody
from sympy.physics.mechanics import Lagrangian, LagrangesMethod
from sympy import symbols

q = q1, q2, q3, q4, q5, q6 = dynamicsymbols('q1:7')
m, g, k, px, Ip = symbols('m g k px Ip')
N = ReferenceFrame('N')
B = N.orientnew('B', 'body', [q4, q5, q6], 'zyx')

A = Point('A')
S = A.locatenew('S', q1*N.x + q2*N.y + q3*N.z)
P = S.locatenew('P', px*B.x)
A.set_vel(N, 0)
S.set_vel(N, S.pos_from(A).dt(N))
P.v2pt_theory(S, N, B)

Ixx = Ip/2
Iyy = Ip/2
Izz = Ip
I = inertia(B, Ixx, Iyy, Izz, 0, 0, 0)
rb = RigidBody('rb', S, B, m, (I, S))
rb.set_potential_energy(
        -m*g*(rb.masscenter.pos_from(A) & N.z) +
        k/2*(P.pos_from(A)).magnitude()**2)

L = Lagrangian(N, rb)
print('{} = {}\n'.format('L', msprint(L)))

lm = LagrangesMethod(L, q)
print(msprint(lm.form_lagranges_equations()))
