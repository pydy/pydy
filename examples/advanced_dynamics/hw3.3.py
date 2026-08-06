#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced Dynamics Homework 3 Exercise 3."""

from sympy.physics.vector import dot, dynamicsymbols
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import msprint
from sympy.physics.mechanics import Point
from sympy import pi, solve, symbols, simplify
from sympy import acos, sin


# define euler angle symbols and reference frames
q1, q2 = dynamicsymbols('q1 q2')
q1d, q2d = dynamicsymbols('q1 q2', 1)
theta, r, R = symbols('θ r R', real=True, positive=True)

theta_val = pi/3
N = ReferenceFrame('N')
#B = N.orientnew('B', 'body', [q1, theta, q2], 'zxz')
F1 = N.orientnew('F1', 'axis', [q1, N.z])
F2 = F1.orientnew('F2', 'axis', [theta, F1.x])
B = F2.orientnew('B', 'axis', [q2, F2.z])

print('rotation matrix from frame N to frame B:')
C = B.dcm(N).subs(theta, theta_val)
print(msprint(C))

# velocity of the disk at the point of contact with the ground is not moving
# since the disk rolls without slipping.
pA = Point('pA') # ball bearing A
pB = pA.locatenew('pB', -R*F1.y) # ball bearing B
pC = pB.locatenew('pC', -r*F2.y) # disk ground contact point

##pO.set_vel(N, 0)
pA.set_vel(N, 0)
pA.set_vel(F1, 0)

pB.set_vel(F1, 0)
pB.set_vel(B, 0)
pB.v2pt_theory(pA, N, F1)

pC.v2pt_theory(pB, N, B)
print('\nvelocity of point C in N, v_C_N, at q1 = 0 = ')
print(pC.vel(N).express(N))

print('\ndisc is rolling without slip so v_C_N = 0')
print('solving for {} with {} = {}'.format(msprint(q2d), theta, theta_val))
soln = solve(pC.vel(N).magnitude(), q2d)
omega2_val = simplify(soln[0].subs(theta, theta_val))
print('{} = {}'.format(msprint(q2d), msprint(soln[0])))
print('{} = {}'.format(msprint(q2d), msprint(omega2_val)))

print('\nthe instantaneous axis of rotation of the disk is simply the')
print('angular velocity vector with respect to frame N')
omega = B.ang_vel_in(N).subs({theta: theta_val, q2d: omega2_val}).express(N)
print('ω = {}'.format(omega))

print('\nto be horizontal, dot(N.z, ω) = 0')
R_val = solve(dot(N.z, omega), R)[0]
print('{} = {} = {}'.format(R, R_val, R_val.subs(theta, theta_val)))

print('\nto find the principal axis of rotation \'a\' and principal '
      'angle of rotation \'φ\'')
print('we can use the equations')
print('    tr(C) = 1 + 2*cos(φ)')
print('    a_x = (c_yz - c_zy)/(2*sin(φ))')
print('    a_y = (c_zx - c_xz)/(2*sin(φ))')
print('    a_z = (c_xy - c_yx)/(2*sin(φ))')
C = C.subs({q1: pi/2, q2:0})
phi = acos((C.trace() - 1)/2)
a_x = (C[1, 2] - C[2, 1])/(2*sin(phi))
a_y = (C[2, 0] - C[0, 2])/(2*sin(phi))
a_z = (C[0, 1] - C[1, 0])/(2*sin(phi))
a = (a_x*N.x + a_y*N.y + a_z*N.z).simplify()
print('\nφ = {} = {}'.format(msprint(phi), phi.n()))
print('a = {}'.format(a))

def eval_vec(vec):
    from sympy.physics.vector import Vector
    v = Vector(0)
    for frame_vectors in vec.args:
        for coeff, uv in zip(list(frame_vectors[0]), frame_vectors[1]):
            v += coeff.n() * uv
    return v

print('  = {}'.format(eval_vec(a)))
