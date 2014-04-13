#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 11.10 from Kane 1985."""

from __future__ import division
from sympy import expand, expand_trig, solve, symbols, trigsimp, sin, cos
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import generalized_active_forces, generalized_inertia_forces
from util import partial_velocities, subs


# Define generalized coordinates, speeds, and constants
q1, q2 = q = dynamicsymbols('q1:3')
q1d, q2d = qd = dynamicsymbols('q1:3', level=1)
u1, u2, u3  = u = dynamicsymbols('u1:4')

mA, mB, mC, t = symbols('mA mB mC t') # masses, time
kA, kB = symbols('kA kB') # radii of gyration
a, a_star, b, b_star, c_star = symbols('a a⁺ b b⁺ c⁺') # distances
#a, a_star, b, b_star, c_star = symbols('a (a*) b (b*) (c*)') # distances
# torque constants
alpha1, alpha2, alpha3 = symbols('α1:4') # exerted on A
beta1, beta2, beta3 = symbols('β1:4') # exerted on C
gamma1, gamma2, gamma3 = symbols('γ1:4') # exerted on A by B
# force constants
P1, P2, P3 = symbols('P1:4') # exerted on A*
Q1, Q2, Q3 = symbols('Q1:4') # exerted on C*
R1, R2, R3 = symbols('R1:4') # exerted on A at P by B
# inertia variables
A1, A3, B1, B3 = symbols('A1 A3 B1 B3')
C1, C2, C3 = symbols('C1:4')

# reference frames
N = ReferenceFrame('N') # also equal to reference frame of cylinder D
A = N.orientnew('A', 'Axis', [q1, N.y]) # counter-weighted crank
B = N.orientnew('B', 'Axis', [q2, -N.y]) # connected rod
C = N.orientnew('C', 'Axis', [0, N.x]) # piston

# points, velocities
pO = Point('O')
pO.set_vel(N, 0)

pA_star = pO.locatenew('A*', -a_star*A.z)
pA_star.v2pt_theory(pO, N, A)

pP = pO.locatenew('P', a*A.z)
pP.v2pt_theory(pO, N, A)

pB_star = pP.locatenew('B*', (b - b_star)*B.z)
pB_star.v2pt_theory(pP, N, B)

# P' is the same point as P but on rigidbody B instead of rigidbody A.
pP_prime = pB_star.locatenew('P\'', -(b - b_star)*B.z)
pP_prime.v2pt_theory(pB_star, N, B)

# Define point Q where B and C intersect.
pQ = pP.locatenew('Q', b*B.z)
pQ.v2pt_theory(pP, N, B)

# Define the distance between points Q, C* as c.
pC_star = pQ.locatenew('C*', c_star*C.z)
pC_star.v2pt_theory(pQ, N, C)

# configuration constraint for q2.
cc = [dot(pC_star.pos_from(pO), N.x)]
cc_map = solve(cc, q2)[1]

# kinematic differential equations
kde = [u1 - q1d, u2 - dot(B.ang_vel_in(N), N.y)]
kde_map = solve(kde, [q1d, q2d])

# velocity constraints
vc = subs([u3 - dot(pC_star.vel(N), N.z), cc[0].diff(t)], kde_map)
vc_map = solve(vc, [u2, u3])

# verify motion constraint equation match text
u2_expected = -a*cos(q1)/(b*cos(q2))*u1
u3_expected = -a/cos(q2)*(sin(q1)*cos(q2) + cos(q1)*sin(q2))*u1
assert trigsimp(vc_map[u2] - u2_expected) == 0
assert trigsimp(vc_map[u3] - u3_expected) == 0

# add the term to get u3 from u1 to kde_map
kde_map[dot(pC_star.vel(N), N.z)] = u3
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

# central inertias, rigid bodies
IA = inertia(A, A1, mA*kA**2, A3)
IB = inertia(B, B1, mB*kB**2, B3)
IC = inertia(C, C1, C2, C3)

# inertia is defined as (central inertia, mass center) for each rigid body
rbA = RigidBody('rbA', pA_star, A, mA, (IA, pA_star))
rbB = RigidBody('rbB', pB_star, B, mB, (IB, pB_star))
rbC = RigidBody('rbC', pC_star, C, mC, (IC, pC_star))
bodies = [rbA, rbB, rbC]

# forces, torques
forces = [(pO, P1*N.x + P2*N.y + P3*N.z),
          (pP, R1*N.x + R2*N.y + R3*N.z),
          (pP_prime, -(R1*N.x + R2*N.y + R3*N.z)),
          (pC_star, Q1*N.x + Q2*N.y + Q3*N.z)]
torques = [(A, alpha1*N.x + alpha2*N.y + alpha3*N.z),
           (A, gamma1*N.x + 0*N.y + gamma3*N.z),
           (B, -(gamma1*N.x + 0*N.y + gamma3*N.z)),
           (C, beta1*N.x + beta2*N.y + beta3*N.z)]

# partial velocities
system = [x for y in bodies for x in [y.masscenter, y.frame]]
system += [f[0] for f in forces + torques]
partials = partial_velocities(system, u, N, kde_map)

# Rewrite the partial velocities of points B*, C*, P'
eq_gen_speed_map = {
        a*u1*cos(q1): -b*u2*cos(q2),
        a*u1*sin(q1): -u3 - a*u1*sin(q2)*cos(q1)/cos(q2),
        (b - b_star)*u2*sin(q2): u3 + a*u1*sin(q1) - b_star*u2*sin(q2)}
for p in [pB_star, pC_star, pP_prime]:
    v = p.vel(N).express(N).subs(kde_map).subs(eq_gen_speed_map)
    partials[p] = dict(zip(u, map(lambda x: v.diff(x, N), u)))

u1d, u2d, u3d = ud = [x.diff(t) for x in u]
for k, v in vc_map.items():
    vc_map[k.diff(t)] = v.diff(t).subs(kde_map).subs(vc_map)

# generalized active/inertia forces
Fr, _ = generalized_active_forces(partials, forces + torques)
Fr_star, _ = generalized_inertia_forces(partials, bodies, kde_map)
print('Fr')
for i, fr in enumerate(Fr, 1):
    print('{0}: {1}'.format(i, msprint(fr)))
print('Fr_star')
for i, fr in enumerate(Fr_star, 1):
    fr = trigsimp(expand(expand_trig(fr)), deep=True, recursive=True)
    print('{0}: {1}'.format(i, msprint(fr)))

# The dynamical equations would be of the form Fr = Fr* (r = 1, 2, 3).
Fr_expected = [
        alpha2 + a*(R1*cos(q1) - R3*sin(q1)),
        b*(R1*cos(q2) + R3*sin(q2)),
        Q3 - R3]
Fr_star_expected = [
        -mA*(a_star**2 + kA**2)*u1d,
        -mB*((b_star**2 + kB**2)*u2d - b_star*sin(q2)*u3d),
        -1*((mB + mC)*u3d - mB*b_star*sin(q2)*u2d + mB*b_star*cos(q2)*u2**2)]
for x, y in zip(Fr, Fr_expected):
    assert expand(x - y) == 0
for x, y in zip(Fr_star, Fr_star_expected):
    assert trigsimp(expand(expand_trig((x - y).subs(vc_map)))) == 0

