#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.12 from Kane 1985."""

from __future__ import division
from sympy import solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot, dynamicsymbols
from util import msprint, subs, partial_velocities, generalized_active_forces


## --- Declare sy3mbols ---
q1 = dynamicsymbols('q1')
u1, u2, u3 = dynamicsymbols('u1:4')
u_prime, R, M, g, e, f, theta = symbols('u\' R, M, g, e, f, theta')
Q1, Q2, Q3 = symbols('Q1, Q2 Q3')

# --- Reference Frames ---
F = ReferenceFrame('F')
P = F.orientnew('P', 'axis', [-theta, F.y])
A = P.orientnew('A', 'axis', [q1, P.x])
A.set_ang_vel(F, u1*A.x + u3*A.z)

## --- define points D, S*, Q on frame A and their velocities ---
pD = Point('D')
pD.set_vel(A, 0)
# u3 will not change v_D_F since wheels are still assumed to roll without slip.
pD.set_vel(F, u2 * A.y)

pS_star = pD.locatenew('S*', e*A.y)
pQ = pD.locatenew('Q', f*A.y - R*A.x)
for p in [pS_star, pQ]:
    p.set_vel(A, 0)
    p.v2pt_theory(pD, F, A)

## --- define partial velocities ---
partials = partial_velocities([pD, pS_star, pQ], [u1, u2, u3],
                              F, express_frame=A)

forces = [(pS_star, -M*g*F.x), (pQ, Q1*A.x + Q2*A.y + Q3*A.z)]
torques = []
Fr, _ = generalized_active_forces(partials, forces + torques)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0} = {1}".format(i, msprint(f)))

F3 = symbols('F3')
fric_Q = Q2*A.y + Q3*A.z
# Q1 is the component of the contact force normal to plane P.
mag_friction_map = {fric_Q.magnitude() : u_prime * Q1}

# friction force points in opposite direction of velocity of Q
vel_Q_F = pQ.vel(F).subs(u3, 0)
eqs = subs([dot(fric_Q.normalize(), A.y) + dot(vel_Q_F.normalize(), A.y),
            dot(fric_Q.normalize(), A.z) + dot(vel_Q_F.normalize(), A.z)],
           mag_friction_map)
Qvals = solve(eqs, [Q2, Q3])

# solve for Q1 in terms of F3 and other variables
Q1_val = solve(F3 - Fr[2].subs(Qvals), Q1)[0]
Qvals[Q1] = Q1_val

print("Contact force components:")
for k in sorted(Qvals.keys()):
    print("{0} = {1}".format(k, msprint(Qvals[k])))
