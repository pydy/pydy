#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 4 Exercise 1."""

from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import inertia, msprint
from sympy.physics.mechanics import Point, RigidBody
from sympy import expand, pi, solve, symbols, simplify
from sympy import acos, sin, cos


q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', 1)
v, m0, m, L, r = symbols('v m0 m L r', real=True, positive=True)

N = ReferenceFrame('N')
C = N.orientnew('C', 'axis', [q1, N.z])
A = C.orientnew('A', 'axis', [q2, C.x])
B = C.orientnew('B', 'axis', [q3, C.x])

pC_star = Point('pC*') # center of mass of rod
pA_star = pC_star.locatenew('pA*', L/2*C.x) # center of disk A
pB_star = pC_star.locatenew('pB*', -L/2*C.x) # center of disk A
pA_hat = pA_star.locatenew('pA^', -r*C.z) # contact point of disk A and ground
pB_hat = pB_star.locatenew('pB^', -r*C.z) # contact point of disk A and ground

pC_star.set_vel(N, v*C.y)
pA_star.v2pt_theory(pC_star, N, C) # pA* and pC* are both fixed in frame C
pB_star.v2pt_theory(pC_star, N, C) # pB* and pC* are both fixed in frame C
pA_hat.v2pt_theory(pA_star, N, A) # pA* and pA^ are both fixed in frame A
pB_hat.v2pt_theory(pB_star, N, B) # pB* and pB^ are both fixed in frame B


I_rod = inertia(C, 0, m0*L**2/12, m0*L**2/12, 0, 0, 0)
rbC = RigidBody('rod_C', pC_star, C, m0, (I_rod, pC_star))

I_discA = inertia(A, m*r**2/2, m*r**2/4, m*r**2/4, 0, 0, 0)
rbA = RigidBody('disc_A', pA_star, A, m, (I_discA, pA_star))

I_discB = inertia(B, m*r**2/2, m*r**2/4, m*r**2/4, 0, 0, 0)
rbB = RigidBody('disc_B', pB_star, B, m, (I_discB, pB_star))

print('omega_A_N = {}'.format(msprint(A.ang_vel_in(N).express(C))))
print('v_pA*_N = {}'.format(msprint(pA_hat.vel(N))))

qd_val = solve([dot(pA_hat.vel(N), C.y), dot(pB_hat.vel(N), C.y)],
               [q2d, q3d])
print(msprint(qd_val))

print('T_A = {}'.format(msprint(simplify(rbA.kinetic_energy(N).subs(qd_val)))))
print('T_B = {}'.format(msprint(simplify(rbB.kinetic_energy(N).subs(qd_val)))))

print('T = {}'.format(msprint(expand(simplify(
    rbA.kinetic_energy(N).subs(qd_val) +
    rbB.kinetic_energy(N).subs(qd_val) +
    rbC.kinetic_energy(N))))))

#T = rb_disc.kinetic_energy(N).subs({theta: theta_val, q2d: q2d_val})
#print('T = {}'.format(msprint(simplify(T))))
#
#t = symbols('t')
#dT = T.diff(symbols('t'))
#print('dT/dt = {} = 0'.format(msprint(simplify(dT))))
