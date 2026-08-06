#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 2 Exercise 4."""

from sympy.physics.vector import ReferenceFrame, Vector
from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.mechanics import mprint, msprint
from sympy import solve

# define euler angle symbols and reference frames
ea = dynamicsymbols('α β γ')
A = ReferenceFrame('A')
A_1 = A.orientnew('A_1', 'axis', [ea[0], A.x])
A_2 = A_1.orientnew('A_2', 'axis', [ea[1], A_1.y])
B = A_2.orientnew('B', 'axis', [ea[2], A_2.z])

# display the rotation matrix C from frame A to frame B
print('Rotation matrix C:')
mprint(B.dcm(A))

# define angular velocity vector
w = dynamicsymbols('ω_x ω_y ω_z')
omega_A = sum((a*uv for a, uv in zip(w, A)), Vector(0))
omega_B = sum((a*uv for a, uv in zip(w, B)), Vector(0))

# display angular velocity vector omega
print('\nangular velocity:')
omega = B.ang_vel_in(A)
print('ω = {}'.format(msprint(omega)))
print('ω_A = {}'.format(msprint(omega.express(A))))
print('ω_B = {}'.format(msprint(omega.express(B))))

# solve for alpha, beta, gamma in terms of angular velocity components
dea = [a.diff() for a in ea]
for omega_fr, frame in [(omega_A, A), (omega_B, B)]:
    eqs = [dot(omega - omega_fr, uv) for uv in frame]
    soln = solve(eqs, dea)
    print('\nif angular velocity components are taken to be in frame '
          '{}'.format(frame))
    for a in dea:
        print('{} = {}'.format(msprint(a), msprint(soln[a])))
