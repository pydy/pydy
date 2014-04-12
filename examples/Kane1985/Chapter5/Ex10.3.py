#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 10.3 from Kane 1985."""

from __future__ import division
from sympy import collect, expand, sin, cos, pi, radsimp, solve, sqrt, symbols
from sympy.physics.mechanics import ReferenceFrame, RigidBody, Point
from sympy.physics.mechanics import dot, dynamicsymbols, inertia, msprint
from util import subs


q1, q2, q3, q4, q5, q6, q7 = q = dynamicsymbols('q1:8')
u1, u2, u3, u4, u5, u6, u7 = u = dynamicsymbols('q1:8', level=1)

M, J, I11, I22, m, r, b = symbols('M J I11 I22 m r b', real=True, positive=True)
omega, t = symbols('ω t')

theta = 30 * pi / 180 # radians
b = r * (1 + sin(theta)) / (cos(theta) - sin(theta))
# Note: using b as found in Ex3.10. Pure rolling between spheres and race R
# is likely a typo and should be between spheres and cone C.

# define reference frames
R = ReferenceFrame('R') # fixed race rf, let R.z point upwards
A = R.orientnew('A', 'axis', [q7, R.z]) # rf that rotates with S* about R.z
# B.x, B.z are parallel with face of cone, B.y is perpendicular
B = A.orientnew('B', 'axis', [-theta, A.x])
S = ReferenceFrame('S')
S.set_ang_vel(A, u1*A.x + u2*A.y + u3*A.z)
C = ReferenceFrame('C')
C.set_ang_vel(A, u4*B.x + u5*B.y + u6*B.z)
#C.set_ang_vel(A, u4*A.x + u5*A.y + u6*A.z)

# define points
pO = Point('O')
pS_star = pO.locatenew('S*', b*A.y)
pS_hat = pS_star.locatenew('S^', -r*B.y) # S^ touches the cone
pS1 = pS_star.locatenew('S1', -r*A.z) # S1 touches horizontal wall of the race
pS2 = pS_star.locatenew('S2', r*A.y) # S2 touches vertical wall of the race

# Since S is rolling against R, v_S1_R = 0, v_S2_R = 0.
pO.set_vel(R, 0)
pS_star.v2pt_theory(pO, R, A)
pS1.v2pt_theory(pS_star, R, S)
pS2.v2pt_theory(pS_star, R, S)
vc = [dot(p.vel(R), basis) for p in [pS1, pS2] for basis in R]

# Since S is rolling against C, v_S^_C = 0.
pO.set_vel(C, 0)
pS_star.v2pt_theory(pO, C, A)
pS_hat.v2pt_theory(pS_star, C, S)
vc += [dot(pS_hat.vel(C), basis) for basis in A]

# Cone has only angular velocity ω in R.z direction.
vc += [dot(C.ang_vel_in(R), basis) for basis in [R.x, R.y]]
vc += [omega - dot(C.ang_vel_in(R), R.z)]

vc_map = solve(vc, u)

# cone rigidbody
I_C = inertia(A, I11, I22, J)
rbC = RigidBody('rbC', pO, C, M, (I_C, pO))
# sphere rigidbody
I_S = inertia(A, 2*m*r**2/5, 2*m*r**2/5, 2*m*r**2/5)
rbS = RigidBody('rbS', pS_star, S, m, (I_S, pS_star))

# kinetic energy
K = radsimp(expand((rbC.kinetic_energy(R) +
                    4*rbS.kinetic_energy(R)).subs(vc_map)))
print('K = {0}'.format(msprint(collect(K, omega**2/2))))

K_expected = (J + 18*m*r**2*(2 + sqrt(3))/5) * omega**2/2
#print('K_expected = {0}'.format(msprint(collect(expand(K_expected),
#                                                omega**2/2))))
assert expand(K - K_expected) == 0
