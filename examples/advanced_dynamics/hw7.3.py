#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 7 Exercise 3."""

from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import msprint
from sympy.physics.mechanics import Point, Particle
from sympy.physics.mechanics import Lagrangian, LagrangesMethod
from sympy import symbols, sin, cos

s, theta, h = dynamicsymbols('s θ h') # s, theta, h
sd, thetad, hd = dynamicsymbols('s θ h', 1)
m, g, l = symbols('m g l')
N = ReferenceFrame('N')

# part a
r1 = s*N.x
r2 = (s + l*cos(theta))*N.x + l*sin(theta)*N.y

O = Point('O')
p1 = O.locatenew('p1', r1)
p2 = O.locatenew('p2', r2)

O.set_vel(N, 0)
p1.set_vel(N, p1.pos_from(O).dt(N))
p2.set_vel(N, p2.pos_from(O).dt(N))

P1 = Particle('P1', p1, 2*m)
P2 = Particle('P2', p2, m)

P1.set_potential_energy(0)
P2.set_potential_energy(P2.mass * g * (p2.pos_from(O) & N.y))

L1 = Lagrangian(N, P1, P2)
print('{} = {}'.format('L1', msprint(L1)))

lm1 = LagrangesMethod(L1, [s, theta])
lm1.form_lagranges_equations()
rhs = lm1.rhs()
t = symbols('t')
print('{} = {}'.format(msprint(sd.diff(t)), msprint(rhs[2].simplify())))
print('{} = {}\n'.format(msprint(thetad.diff(t)), msprint(rhs[3].simplify())))

# part b
r1 = s*N.x + h*N.y
r2 = (s + l*cos(theta))*N.x + (h + l*sin(theta))*N.y

p1 = O.locatenew('p1', r1)
p2 = O.locatenew('p2', r2)
p1.set_vel(N, p1.pos_from(O).dt(N))
p2.set_vel(N, p2.pos_from(O).dt(N))

P1 = Particle('P1', p1, 2*m)
P2 = Particle('P2', p2, m)

P1.set_potential_energy(P1.mass * g * (p1.pos_from(O) & N.y))
P2.set_potential_energy(P2.mass * g * (p2.pos_from(O) & N.y))

L2 = Lagrangian(N, P1, P2)
print('{} = {}'.format('L2', msprint(L2)))

lm2 = LagrangesMethod(L2, [s, theta, h], hol_coneqs = [h])
lm2.form_lagranges_equations()
rhs = lm2.rhs()
print('{} = {}'.format(msprint(sd.diff(t)), msprint(rhs[3].simplify())))
print('{} = {}'.format(msprint(thetad.diff(t)), msprint(rhs[4].simplify())))
print('{} = {}'.format(msprint(hd.diff(t)), msprint(rhs[5].simplify())))
print('{} = {}'.format('λ', msprint(rhs[6].simplify())))
