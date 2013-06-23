#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.18 from Kane 1985.
"""

from __future__ import division
from sympy import symbols
from sympy.physics.mechanics import ReferenceFrame
from sympy.physics.mechanics import cross, dot, dynamicsymbols, inertia
from util import msprint

print("\n part a")
Ia, Ib, Ic, Iab, Ibc, Ica, t = symbols('Ia Ib Ic Iab Ibc Ica t')
omega = dynamicsymbols('omega')
N = ReferenceFrame('N')

# I = (I11 * N.x + I12 * N.y + I13 * N.z) N.x + 
#     (I21 * N.x + I22 * N.y + I23 * N.z) N.y + 
#     (I31 * N.x + I32 * N.y + I33 * N.z) N.z

# definition of T* is:
# T* = -dot(alpha, I) - dot(cross(omega, I), omega)
ang_vel = omega * N.x
I = inertia(N, Ia, Ib, Ic, Iab, Ibc, Ica)

T_star = -dot(ang_vel.diff(t, N), I) - dot(cross(ang_vel, I), ang_vel)
print(msprint(T_star))

print("\n part b")
I11, I22, I33, I12, I23, I31 = symbols('I11 I22 I33 I12 I23 I31')
omega1, omega2, omega3 = dynamicsymbols('omega1:4')
B = ReferenceFrame('B')
ang_vel = omega1 * B.x + omega2 * B.y + omega3 * B.z
I = inertia(B, I11, I22, I33, I12, I23, I31)
T_star = -dot(ang_vel.diff(t, B), I) - dot(cross(ang_vel, I), ang_vel)
print(msprint(T_star))

