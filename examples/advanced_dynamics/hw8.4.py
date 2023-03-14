#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 7 Exercise 3."""

from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import msprint, inertia
from sympy.physics.mechanics import RigidBody, Point
from sympy.physics.mechanics import Lagrangian, LagrangesMethod
from sympy import symbols, sin, cos

q = q1, q2, q3 = dynamicsymbols('q1:4') # x, y, theta
qd = q1d, q2d, q3d = dynamicsymbols('q1:4', 1)
t, g, m, l, w, f, v0 = symbols('t g m l w f v0')
Fx, Fy = symbols('Fx Fy')
values = {
    g: 9.81,
    m: 20,
    l: 2,
    w: 1,
    f: 2,
    v0: 20}

N = ReferenceFrame('N')
B = N.orientnew('B', 'axis', [q3, N.z])


O = Point('O')
S = O.locatenew('S', q1*N.x + q2*N.y)
S.set_vel(N, S.pos_from(O).dt(N))

#Is = m/12*(l**2 + w**2)
Is = symbols('Is')
I = inertia(B, 0, 0, Is, 0, 0, 0)
rb = RigidBody('rb', S, B, m, (I, S))
rb.set_potential_energy(0)


L = Lagrangian(N, rb)
lm = LagrangesMethod(
    L, q, nonhol_coneqs = [q1d*sin(q3) - q2d*cos(q3) + l/2*q3d])
lm.form_lagranges_equations()
rhs = lm.rhs()
print('{} = {}'.format(msprint(q1d.diff(t)),
    msprint(rhs[3].simplify())))
print('{} = {}'.format(msprint(q2d.diff(t)),
    msprint(rhs[4].simplify())))
print('{} = {}'.format(msprint(q3d.diff(t)),
    msprint(rhs[5].simplify())))
print('{} = {}'.format('λ', msprint(rhs[6].simplify())))

