#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.13 from Kane 1985."""

from __future__ import division
from sympy import solve, symbols, trigsimp
from sympy.physics.mechanics import Point, ReferenceFrame, RigidBody
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import inertia_coefficient_matrix, partial_velocities


q1, q2, q3 = q = dynamicsymbols('q1:4')
q1d, q2d, q3d = qd = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = u = dynamicsymbols('u1:4')
m, IB11, IB22, IB33 = symbols('m IB11 IB22 IB33')

# reference frames
A = ReferenceFrame('A')
B = A.orientnew('B', 'body', [q1, q2, q3], 'xyz')

# define a point fixed in both A and B
pO = Point('O')
pO.set_vel(A, 0)

# define the rigid body B
I = inertia(B, IB11, IB22, IB33)
rb = RigidBody('rbB', pO, B, m, (I, pO))


def coupled_speeds(ic_matrix, partials):
    ulist = partials.ulist
    print('dynamically coupled speeds:')
    found = False
    for i, r in enumerate(ulist):
        for j, s in enumerate(ulist[i + 1:], i + 1):
            if trigsimp(ic_matrix[i, j]) != 0:
                found = True
                print('{0} and {1}'.format(msprint(r), msprint(s)))
    if not found:
        print('None')


def find_coupled_speeds(kde_map):
    partials = partial_velocities([rb.frame, rb.masscenter], [u1, u2, u3],
                                  A, kde_map)
    M = trigsimp(inertia_coefficient_matrix([rb], partials))
    coupled_speeds(M, partials)
    return M

print('part a')
kde_map = dict(zip(qd, u))
M = find_coupled_speeds(kde_map)
assert M[0, 1] != 0
assert M[0, 2] != 0
assert M[1, 2] == 0

print('\npart b')
kde = [x - y for x, y in zip(u, [dot(B.ang_vel_in(A), basis) for basis in B])]
kde_map = solve(kde, qd)
M = find_coupled_speeds(kde_map)
assert M[0, 1] == 0
assert M[0, 2] == 0
assert M[1, 2] == 0
