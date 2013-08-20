#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.3 from Kane 1985."""

from __future__ import division
from sympy import collect, expand, sin, cos, pi, solve, symbols, trigsimp, sqrt, expand_trig
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import subs


q1, q2, q3, q4 = q = dynamicsymbols('q1:5')
q1d, q2d, q3d, q4d = qd = dynamicsymbols('q1:5', level=1)
u1, u2, u3, u4 = u = dynamicsymbols('u1:5')

#theta, omega, t = symbols('θ ω t')
omega, t = symbols('ω t')
#theta, omega, m, r, t = symbols('theta omega m r t') # real=True, positive=True)
# M, I11, I33 should not show up kinetic energy
M, J, I11, I33, m, r, theta = symbols('M J I11 I33 m r θ', real=True, positive=True)
b = symbols('b', real=True, positive=True)

#theta = 30 * pi/180 # radians (30 degrees)

## --- Define ReferenceFrames ---
R = ReferenceFrame('R') # fixed race rf
# R.y is parallel with shaft axis, R.x point outwards from shaft axis to race
C = R.orientnew('C', 'axis', [omega * t, R.y])
#S_prime = C.orientnew('S_prime', 'axis', [q4, R.y])
S_prime = C.orientnew('S\'', 'axis', [q4, C.y])
S = S_prime.orientnew('S', 'body', [q1, q2, q3], 'xyz')
#S = R.orientnew('S', 'body', [q1, q2, q3], 'xyz')

## --- Define Points and their velocities ---
#b = r*(1 + sin(theta))/(cos(theta) - sin(theta))

pO = Point('O') # point where the shaft/cone touches the race
pO.set_vel(R, 0)
pO.set_vel(C, 0)
pO.set_vel(S_prime, 0)

# center of a sphere S
pS_star = pO.locatenew('S*', b*S_prime.x + r*S_prime.y)
pS_star.set_vel(S_prime, 0)
pS_star.set_vel(S, 0)
#pS0 = pS_star.locatenew('S0', r*(-sin(theta)*S.x +
#                                 cos(theta)*S.y))
## point S1 is on S touching vertical face of R
#pS1 = pS_star.locatenew('S1', r*S.x)
## point S2 is on S touching horizontal face of R
#pS2 = pS_star.locatenew('S1', -r*S.y)
pS0 = pS_star.locatenew('S0', r*(-cos(theta)*S_prime.x +
                                 sin(theta)*S_prime.y))
pS0.set_vel(S, 0)
# point S1 is on S touching vertical face of R
pS1 = pS_star.locatenew('S1', r*S_prime.x)
pS1.set_vel(S, 0)
# point S2 is on S touching horizontal face of R
pS2 = pS_star.locatenew('S2', -r*S_prime.y)
pS2.set_vel(S, 0)

# calculate velocities
pS_star.v2pt_theory(pO, C, S_prime)
pS_star.v2pt_theory(pO, R, S_prime)
for rf in [S_prime, C, R]:
    for p in [pS0, pS1, pS2]:
        p.v2pt_theory(pS_star, rf, S)

## --- Expressions for generalized speeds u1, u2, u3, u4, u5 ---
kde_map = dict(zip(qd, u))
#vc = [pS0.vel(C), pS1.vel(R), pS2.vel(R)] # rolling constraints
# pS0.vel(C) = 0: used to solve for b
vc = [dot(v, b) for v in [pS1.vel(R), pS2.vel(R)] for b in R]
# pure rolling
vc += [dot(S.ang_vel_in(R), R.x) + dot(S.ang_vel_in(R), R.y)]
vc = map(lambda x: trigsimp(expand_trig(x)), vc)
for cons in vc:
    print('\n{0}'.format(msprint(cons)))
vc_map = solve(vc, qd)
print(msprint(vc_map))
print(msprint(simplify(S.ang_vel_in(C).express(C))))
if not vc_map:
    import sys
    sys.exit(1)

# cone rigidbody
I_C = inertia(C, I11, J, I33)
rbC = RigidBody('rbC', pO, C, M, (I_C, pO))
# sphere rigidbody
I_S = inertia(S, 2*m*r**2/3, 2*m*r**2/3, 2*m*r**2/3)
rbS = RigidBody('rbS', pS_star, S, m, (I_S, pS_star))

# kinetic energy
K = trigsimp((rbC.kinetic_energy(R) + 4*rbS.kinetic_energy(R))).subs(vc_map)
print('K = {0}'.format(msprint(collect(K, omega**2/2))))

K_expected = (J + 18*m*r**2*(2 + sqrt(3))/5) * omega**2/2
assert expand(K - K_expected) == 0
