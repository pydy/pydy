#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 5 Exercise 2."""

from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import inertia, msprint
from sympy.physics.mechanics import Point, RigidBody
from sympy import pi, solve, symbols, simplify
from sympy import acos, sin, cos


# 2a
q1, q2 = dynamicsymbols('q1 q2')
q1d, q2d = dynamicsymbols('q1 q2', 1)
r, R, m, g = symbols('r R m g', real=True, positive=True)
theta = pi/3
q2d_val = (-(R + r*cos(theta))/r*q1d)
vals = {q2d: q2d_val}

N = ReferenceFrame('N')
F1 = N.orientnew('F1', 'axis', [q1, N.z])
F2 = F1.orientnew('F2', 'axis', [theta, F1.x])
B = F2.orientnew('B', 'axis', [q2, F2.z])

# bearing A
pA = Point('A')
pA.set_vel(N, 0)
pA.set_vel(F1, 0)

# bearing B, center of mass of disc
pB = pA.locatenew('pB', -R*F1.y)
pB.set_vel(B, 0)
pB.v2pt_theory(pA, N, F1)
print('v_B_N = {}'.format(msprint(pB.vel(N))))

Ixx = m*r**2/4
Iyy = m*r**2/4
Izz = m*r**2/2
I_disc = inertia(F2, Ixx, Iyy, Izz, 0, 0, 0)
rb_disc = RigidBody('disc', pB, B, m, (I_disc, pB))
H = rb_disc.angular_momentum(pB, N).subs(vals).express(F2).simplify()
print("H about B in frame N = {}".format(msprint(H)))

#2b
# disc/ground contact point
pC = pB.locatenew('pC', -r*F2.y)

fAx, fAy, fAz, fBx, fBy, fBz = symbols('fAx fAy fAz fBx fBy fBz')
fCx, fCy, fCz = symbols('fCx fCy fCz')
mAx, mAy, mAz, mBx, mBy, mBz = symbols('mAx mAy mAz mBx mBy mBz')

# forces on rod, disc
fA = fAx*F1.x + fAy*F1.y + fAz*F1.z # force exerted on rod at point A
fB = fBx*F2.x + fBy*F2.y + fBz*F2.z # force exerted on disc by rod at point B
fC = fCx*F2.x + fCy*F2.y + fCz*F2.z # force exerted on ground by disc at point C
mA = mAx*F1.x + mAy*F1.y + mAz*F1.z # moment exerted on rod A from axis
mB = mBx*F2.x + mBy*F2.y + mBz*F2.z # moment exerted on disc by rod
constraints = [
        mA & F1.z, # no moment due to bearing A
        fA & F1.z, # no force due to bearing A
        mB & F2.z, # no moment due to bearing B
        fB & F2.z, # no force due to bearing B
        (fA - fB) & F1.x, # rod has no mass
        (fA - fB) & F1.y,
        (fA - fB) & F1.z,
        (mA - mB - (pB.pos_from(pA) ^ fB)) & F1.x, # rod has no inertia
        (mA - mB - (pB.pos_from(pA) ^ fB)) & F1.y,
        (mA - mB - (pB.pos_from(pA) ^ fB)) & F1.z,
        ]

print('\nconstraints due to bearings and rod')
for i, c in enumerate(constraints):
    print('{} {}'.format(i, msprint(c)))
soln = solve(constraints,
        [fAx, fAy, fAz, fBx, fBy, fBz, mAx, mAy, mAz, mBx, mBy, mBz])
print('\nmoments and forces after constraints applied:')
print(soln)

print('\nfA = {}'.format(msprint(fA.subs(soln))))
print('fB = {}'.format(msprint(fB.subs(soln))))
print('fC = {}'.format(msprint(fC.subs(soln))))
print('mA = {}'.format(msprint(mA.subs(soln))))
print('mB = {}'.format(msprint(mB.subs(soln))))

#2c
a_B_N = pB.acc(N).express(F1)
print('\na_B_N = {}'.format(msprint(a_B_N)))
dH = H.dt(N)
print('dH/dt = {}'.format(dH))

forces_disc = fB - fC - m*g*N.z
moments_disc = mB - (pC.pos_from(pB) ^ fC)
constraints = [
        ((forces_disc - m*a_B_N) & F1.x).subs(soln),
        ((forces_disc - m*a_B_N) & F1.y).subs(soln),
        ((forces_disc - m*a_B_N) & F1.z).subs(soln),
        ((moments_disc - dH) & F2.x).subs(soln),
        ((moments_disc - dH) & F2.y).subs(soln),
        #((moments_disc - dH) & F2.z).subs(soln), # system overconstrained
        ]
print('constraints')
for i, c in enumerate(constraints):
    print('{} {}'.format(i, msprint(c)))

soln2 = solve(constraints, [mAx, mBy, fCx, fCy, fCz])
print(msprint(soln2))

print('\nforce and moment components:')
for fm in [fA, mA]:
    for uv in F1:
        component = fm & uv
        print('{} = {}'.format(component,
            msprint(component.subs(soln).subs(soln2))))

for fm in [fB, fC, mB]:
    for uv in F2:
        component = fm & uv
        print('{} = {}'.format(component,
            msprint(component.subs(soln).subs(soln2))))
