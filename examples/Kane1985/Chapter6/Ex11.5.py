#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 11.5 from Kane 1985."""

from __future__ import division
from sympy import expand, solve, symbols, trigsimp
from sympy import sin, tan, pi
from sympy.physics.mechanics import Point, ReferenceFrame, RigidBody
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import generalized_active_forces, generalized_inertia_forces
from util import partial_velocities, subs


g, m, Px, Py, Pz, R, t = symbols('g m Px Py Pz R t')
q1, q2, q3, q4, q5 = q = dynamicsymbols('q1:6')
qd = dynamicsymbols('q1:6', level=1)
u1, u2, u3, u4, u5 = u = dynamicsymbols('u1:6')

# reference frames
A = ReferenceFrame('A')
B_prime = A.orientnew('B_prime', 'Axis', [q1, A.z])
B = B_prime.orientnew('B', 'Axis', [pi/2 - q2, B_prime.x])
C = B.orientnew('C', 'Axis', [q3, B.z])

# points, velocities
pO = Point('O')
pO.set_vel(A, 0)

# R is the point in plane H that comes into contact with disk C.
pR = pO.locatenew('R', q4*A.x + q5*A.y)
pR.set_vel(A, pR.pos_from(pO).dt(A))
pR.set_vel(B, 0)

# C^ is the point in disk C that comes into contact with plane H.
pC_hat = pR.locatenew('C^', 0)
pC_hat.set_vel(C, 0)

# C* is the point at the center of disk C.
pC_star = pC_hat.locatenew('C*', R*B.y)
pC_star.set_vel(C, 0)
pC_star.set_vel(B, 0)

# calculate velocities in A
pC_star.v2pt_theory(pR, A, B)
pC_hat.v2pt_theory(pC_star, A, C)

# kinematic differential equations
kde = [x - y for x, y in zip(
        [dot(C.ang_vel_in(A), basis) for basis in B] + qd[3:],
        u)]
kde_map = solve(kde, qd)
# include second derivatives in kde map
for k, v in kde_map.items():
    kde_map[k.diff(t)] = v.diff(t)

vc = map(lambda x: dot(pC_hat.vel(A), x), [A.x, A.y])
vc_map = solve(subs(vc, kde_map), [u4, u5])

# define disc rigidbody
I_C = inertia(C, m*R**2/4, m*R**2/4, m*R**2/2)
rbC = RigidBody('rbC', pC_star, C, m, (I_C, pC_star))

# forces
R_C_hat = Px*A.x + Py*A.y + Pz*A.z
R_C_star = -m*g*A.z
forces = [(pC_hat, R_C_hat), (pC_star, R_C_star)]

# partial velocities
bodies = [rbC]
system = ([i.masscenter for i in bodies] + [i.frame for i in bodies] +
          list(zip(*forces)[0]))
partials = partial_velocities(system, [u1, u2, u3], A, kde_map, vc_map)

# generalized active forces
Fr, _ = generalized_active_forces(partials, forces)
Fr_star, _ = generalized_inertia_forces(partials, bodies, kde_map, vc_map)

# dynamical equations
dyn_eq = subs([x + y for x, y in zip(Fr, Fr_star)], kde_map)
u1d, u2d, u3d = ud = [x.diff(t) for x in [u1, u2, u3]]
dyn_eq_map = solve(dyn_eq, ud)

for x in ud:
    print('{0} = {1}'.format(msprint(x),
                             msprint(trigsimp(dyn_eq_map[x]))))

u1d_expected = (u2**2*tan(q2) - 6*u2*u3 -4*g*sin(q2)/R)/5
u2d_expected = 2*u3*u1 - u1*u2*tan(q2)
u3d_expected = 2*u1*u2/3
assert trigsimp(expand(dyn_eq_map[u1d] - u1d_expected)) == 0
assert trigsimp(expand(dyn_eq_map[u2d] - u2d_expected)) == 0
assert trigsimp(expand(dyn_eq_map[u3d] - u3d_expected)) == 0
