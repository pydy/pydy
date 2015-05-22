#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 6.6 from Kane 1985."""

from __future__ import division
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia, inertia_of_point_mass
from sympy.physics.mechanics import dot
from sympy import symbols
from sympy import S

m = symbols('m')
m_val = 12

N = ReferenceFrame('N')
pO = Point('O')
pBs = pO.locatenew('B*', -3*N.x + 2*N.y - 4*N.z)

I_B_O = inertia(N, 260*m/m_val, 325*m/m_val, 169*m/m_val,
                72*m/m_val, 96*m/m_val, -144*m/m_val)
print("I_B_rel_O = {0}".format(I_B_O))

I_Bs_O = inertia_of_point_mass(m, pBs.pos_from(pO), N)
print("\nI_B*_rel_O = {0}".format(I_Bs_O))

I_B_Bs = I_B_O - I_Bs_O
print("\nI_B_rel_B* = {0}".format(I_B_Bs))

pQ = pO.locatenew('Q', -4*N.z)
I_Bs_Q = inertia_of_point_mass(m, pBs.pos_from(pQ), N)
print("\nI_B*_rel_Q = {0}".format(I_Bs_Q))

I_B_Q = I_B_Bs + I_Bs_Q
print("\nI_B_rel_Q = {0}".format(I_B_Q))

# n_a is a vector parallel to line PQ
n_a = S(3)/5 * N.x - S(4)/5 * N.z
I_a_a_B_Q = dot(dot(n_a, I_B_Q), n_a)
print("\nn_a = {0}".format(n_a))
print("\nI_a_a_B_Q = {0} = {1}".format(I_a_a_B_Q, I_a_a_B_Q.subs(m, m_val)))
