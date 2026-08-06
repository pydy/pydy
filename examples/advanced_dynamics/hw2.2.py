#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 2 Exercise 2."""

from sympy.physics.vector import ReferenceFrame, Vector
from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.mechanics import mprint, msprint
from sympy import pi

# define euler angle symbols and reference frames
ea = dynamicsymbols('α β γ') # zyx
s_i = dynamicsymbols('s_x s_y s_z')
p_i = dynamicsymbols('p_x p_y p_z')
A = ReferenceFrame('A')
B = A.orientnew('B', 'body', ea, 'zyx')

# display the rotation matrix C from frame A to frame B
print('Rotation matrix C:')
mprint(B.dcm(A))

vals = [pi, pi/3, -pi/4]
print('\nsubstitute in euler angle values')
print('    α = {}'.format(vals[0]))
print('    β = {}'.format(vals[1]))
print('    γ = {}'.format(vals[2]))
ea_dict = dict(zip(ea, vals))
R = B.dcm(A).subs(ea_dict)
print('R = ')
mprint(R)

p = s_i[0]*A.x + s_i[1]*A.y + s_i[2]*A.z
p += p_i[0]*B.x + p_i[2]*B.z
print('\nvector p = ')
mprint(p.express(A).subs(ea_dict))

print('\nmagnitude of p = ')
mprint(p.express(A).subs(ea_dict).magnitude())
