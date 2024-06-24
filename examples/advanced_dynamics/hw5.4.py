#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 5 Exercise 4."""

from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import inertia, msprint
from sympy.physics.mechanics import Point, RigidBody
from sympy import pi, solve, symbols, simplify
from sympy import acos, sin, cos


# define euler angle symbols and reference frames
q1, q2 = dynamicsymbols('q1 q2')
q1d, q2d = dynamicsymbols('q1 q2', 1)
q1dd, q2dd = dynamicsymbols('q1 q2', 2)
m, Ia, It = symbols('m Ia It', real=True, positive=True)

N = ReferenceFrame('N')
F = N.orientnew('F', 'axis', [q1, N.z]) # gimbal frame
B = F.orientnew('B', 'axis', [q2, F.x]) # flywheel frame

P = Point('P')
P.set_vel(N, 0)
P.set_vel(F, 0)
P.set_vel(B, 0)

I = inertia(F, Ia, It, It, 0, 0, 0)
rb = RigidBody('flywheel', P, B, m, (I, P))

H = rb.angular_momentum(P, N)
print('H_P_N = {}\n    = {}'.format(H.express(N), H.express(F)))
dH = H.dt(N)
print('d^N(H_P_N)/dt = {}'.format(dH.express(F).subs(q2dd, 0)))

print('\ndH/dt = M')
print('M = {}'.format(dH.express(F).subs(q2dd, 0)))

print('\ncalculation using euler angles')
t = symbols('t')
omega = F.ang_vel_in(N)
wx = omega & F.x
wy = omega & F.y
wz = omega & F.z
s = B.ang_vel_in(F) & F.x
Mx = Ia*(wx + s).diff(t)
My = It*(wy).diff(t) - (It - Ia)*wz*wx + Ia*s*wz
Mz = It*(wz).diff(t) - (It - Ia)*wx*wy + Ia*s*wy
M = Mx*F.x + My*F.y + Mz*F.z
print('M = {}'.format(M.subs(q2dd, 0)))

Fr, a = symbols('Fr a')
f = Fr*F.x
soln = solve([(((2*a*N.z)^f) - M) & F.y], Fr)

print('\nFr = {}'.format(msprint(soln[Fr])))

