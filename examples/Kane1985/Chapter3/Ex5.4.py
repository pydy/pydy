#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 5.4 from Kane 1985."""

from __future__ import division
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot
from sympy import symbols, integrate
from sympy import sin, cos, pi


a, b, R = symbols('a b R', real=True, nonnegative=True)
theta, x, y, z = symbols('theta x y z', real=True)
ab_val = {a: 0.3, b: 0.3}

# R is the radius of the circle, theta is the half angle.
centroid_sector = 2*R*sin(theta) / (3 * theta)

# common R, theta values
theta_pi_4 = {theta: pi/4, R: a}
R_theta_val = {theta: pi/4 * (1 - z/a), R: a}

N = ReferenceFrame('N')
def eval_vec(v):
    vs = v.subs(ab_val)
    return sum(dot(vs, n).evalf() * n for n in N)

# For each part A, B, C, D, define an origin for that part such that the
# centers of mass of each part of the component have positive N.x, N.y,
# and N.z values.
## FIND CENTER OF MASS OF A
vol_A_1 = pi * a**2 * b / 4
vol_A_2 = pi * a**2 * a / 4 / 2
vol_A = vol_A_1 + vol_A_2
pA_O = Point('A_O')
pAs_1 = pA_O.locatenew(
        'A*_1', (b/2 * N.z +
                 centroid_sector.subs(theta_pi_4) * sin(pi/4) * (N.x + N.y)))
pAs_2 = pA_O.locatenew(
        'A*_2', (b * N.z +
                 N.z * integrate((theta*R**2*(z)).subs(R_theta_val),
                                 (z, 0, a)) / vol_A_2 +
                 N.x * integrate((theta*R**2 * cos(theta) *
                                  centroid_sector).subs(R_theta_val),
                                 (z, 0, a)) / vol_A_2 +
                 N.y * integrate((2*R**3/3 * 4*a/pi *
                                  sin(theta)**2).subs(R, a),
                                 (theta, 0, pi/4)) / vol_A_2))
pAs = pA_O.locatenew('A*', ((pAs_1.pos_from(pA_O) * vol_A_1 +
                             pAs_2.pos_from(pA_O) * vol_A_2) /
                            vol_A))
print('A* = {0}'.format(pAs.pos_from(pA_O)))
print('A* = {0}'.format(eval_vec(pAs.pos_from(pA_O))))

## FIND CENTER OF MASS OF B
vol_B_1 = pi*a**2/2
vol_B_2 = a**2 / 2
vol_B  = vol_B_1 + vol_B_2
pB_O = Point('B_O')
pBs_1 = pB_O.locatenew(
        'B*_1', (a*(N.x + N.z) + a/2*N.y +
                 (-N.x + N.z) * (R*sin(theta)/theta *
                                 sin(pi/4)).subs(theta_pi_4)))
pBs_2 = pB_O.locatenew('B*_2', (a*N.y + a*N.z -
                                (a/3 * N.y + a/3 * N.z)))
pBs = pB_O.locatenew('B*', ((pBs_1.pos_from(pB_O) * vol_B_1 +
                             pBs_2.pos_from(pB_O) * vol_B_2) /
                            vol_B))
print('\nB* = {0}'.format(pBs.pos_from(pB_O)))
print('B* = {0}'.format(eval_vec(pBs.pos_from(pB_O))))

## FIND CENTER OF MASS OF C
vol_C_1 = 2 * a**2 * b
vol_C_2 = a**3 / 2
vol_C_3 = a**3
vol_C_4 = -pi*a**3/4
vol_C = vol_C_1 + vol_C_2 + vol_C_3 + vol_C_4
pC_O = Point('C_O')
pCs_1 = pC_O.locatenew('C*_1', (a*N.x + a/2*N.y + b/2*N.z))
pCs_2 = pC_O.locatenew('C*_2', (a*N.x + b*N.z +
                                (a/3*N.x + a/2*N.y + a/3*N.z)))
pCs_3 = pC_O.locatenew('C*_3', (b*N.z + a/2*(N.x + N.y + N.z)))
pCs_4 = pC_O.locatenew(
        'C*_4', ((a + b)*N.z + a/2*N.y +
                 (N.x - N.z)*(centroid_sector.subs(
                         theta_pi_4)*sin(pi/4))))
pCs = pC_O.locatenew('C*', ((pCs_1.pos_from(pC_O)*vol_C_1 +
                             pCs_2.pos_from(pC_O)*vol_C_2 +
                             pCs_3.pos_from(pC_O)*vol_C_3 +
                             pCs_4.pos_from(pC_O)*vol_C_4) /
                            vol_C))
print('\nC* = {0}'.format(pCs.pos_from(pC_O)))
print('C* = {0}'.format(eval_vec(pCs.pos_from(pC_O))))

## FIND CENTER OF MASS OF D
vol_D = pi*a**3/4
pD_O = Point('D_O')
pDs = pD_O.locatenew('D*', (a*N.z + a/2*N.y +
                            (N.x - N.z)*(centroid_sector.subs(
                                     theta_pi_4) * sin(pi/4))))
print('\nD* = {0}'.format(pDs.pos_from(pD_O)))
print('D* = {0}'.format(eval_vec(pDs.pos_from(pD_O))))

## FIND CENTER OF MASS OF ASSEMBLY
pO = Point('O')
pA_O.set_pos(pO, 2*a*N.x - (a+b)*N.z)
pB_O.set_pos(pO, 2*a*N.x - a*N.z)
pC_O.set_pos(pO, -(a+b)*N.z)
pD_O.set_pos(pO, -a*N.z)

density_A = 7800
density_B = 17.00
density_C = 2700
density_D = 8400
mass_A = vol_A * density_A
mass_B = vol_B * density_B
mass_C = vol_C * density_C
mass_D = vol_D * density_D

pms = pO.locatenew('m*', ((pAs.pos_from(pO)*mass_A + pBs.pos_from(pO)*mass_B +
                           pCs.pos_from(pO)*mass_C + pDs.pos_from(pO)*mass_D) /
                          (mass_A + mass_B + mass_C + mass_D)))
print('\nm* = {0}'.format(eval_vec(pms.pos_from(pO))))
