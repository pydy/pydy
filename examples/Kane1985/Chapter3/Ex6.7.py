#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 6.7 from Kane 1985."""

from __future__ import division
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia, inertia_of_point_mass
from sympy.physics.mechanics import cross, dot
from sympy import solve, sqrt, symbols, integrate, diff
from sympy import sin, cos, tan, pi, S, Matrix
from sympy import simplify, Abs

b, m, k_a = symbols('b m k_a', real=True, nonnegative=True)
theta = symbols('theta', real=True)

N = ReferenceFrame('N')
N2 = N.orientnew('N2', 'Axis', [theta, N.z])
pO = Point('O')
pP1s = pO.locatenew('P1*', b/2 * (N.x + N.y))
pP2s = pO.locatenew('P2*', b/2 * (2*N.x + N.y + N.z))
pP3s = pO.locatenew('P3*', b/2 * ((2 + sin(theta))*N.x +
                                  (2 - cos(theta))*N.y +
                                  N.z))

I_1s_O = inertia_of_point_mass(m, pP1s.pos_from(pO), N)
I_2s_O = inertia_of_point_mass(m, pP2s.pos_from(pO), N)
I_3s_O = inertia_of_point_mass(m, pP3s.pos_from(pO), N)
print("\nI_1s_rel_O = {0}".format(I_1s_O))
print("\nI_2s_rel_O = {0}".format(I_2s_O))
print("\nI_3s_rel_O = {0}".format(I_3s_O))


I_1_1s = inertia(N, m*b**2/12, m*b**2/12, 2*m*b**2/12)
I_2_2s = inertia(N, 2*m*b**2/12, m*b**2/12, m*b**2/12)

I_3_3s_N2 = Matrix([[2, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]])
I_3_3s_N = m*b**2/12 * simplify(N.dcm(N2) * I_3_3s_N2 * N2.dcm(N))
I_3_3s = inertia(N, I_3_3s_N[0, 0], I_3_3s_N[1, 1], I_3_3s_N[2, 2],
                    I_3_3s_N[0, 1], I_3_3s_N[1, 2], I_3_3s_N[2, 0])

I_1_O = I_1_1s + I_1s_O
I_2_O = I_2_2s + I_2s_O
I_3_O = I_3_3s + I_3s_O
print("\nI_1_rel_O = {0}".format(I_1_O))
print("\nI_2_rel_O = {0}".format(I_2_O))
print("\nI_3_rel_O = {0}".format(I_3_O))

# assembly inertia tensor is the sum of the inertia tensor of each component
I_B_O = I_1_O + I_2_O + I_3_O
print("\nI_B_rel_O = {0}".format(I_B_O))

# n_a is parallel to line L
n_a = sqrt(2) / 2 * (N.x + N.z)
print("\nn_a = {0}".format(n_a))

# calculate moment of inertia of for point O of assembly about line L
I_a_a_B_O = simplify(dot(n_a, dot(I_B_O, n_a)))
print("\nI_a_a_B_rel_O = {0}".format(I_a_a_B_O))

# use the second value since k_a is non-negative
k_a_val = solve(I_a_a_B_O - 3 * m * k_a**2, k_a)[1]
print("\nk_a = {0}".format(k_a_val))

dk_a_dtheta = diff(k_a_val, theta)
print("\ndk_a/dtheta = {0}".format(dk_a_dtheta))

# solve for theta = 0 using a simplified expression or
# else no solution will be found
soln = solve(3*cos(theta) + 12*sin(theta) - 4*sin(theta)*cos(theta), theta)
# ignore complex values of theta
theta_vals = [s for s in soln if s.is_real]

print("")
for th in theta_vals:
    print("k_a({0}) = {1}".format((th * 180 / pi).evalf(),
                                  k_a_val.subs(theta, th).evalf()))
