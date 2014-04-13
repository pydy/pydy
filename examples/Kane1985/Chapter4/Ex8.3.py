#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.3 from Kane 1985."""

from __future__ import division
#import os
#import sys
#sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
#                             "../../../sympy/sympy/physics/mechanics"))

from sympy import pi, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point, dot, dynamicsymbols
#from sympy.physics.mechanics import KanesMethod
#from kane import KanesMethod
from util import msprint, subs, partial_velocities, generalized_active_forces


g, m, Px, Py, Pz, R, t = symbols('g m Px Py Pz R t')
q = dynamicsymbols('q1:6')
qd = dynamicsymbols('q1:6', level=1)
u = dynamicsymbols('u1:6')

## --- Define ReferenceFrames ---
A = ReferenceFrame('A')
B_prime = A.orientnew('B_prime', 'Axis', [q[0], A.z])
B = B_prime.orientnew('B', 'Axis', [pi/2 - q[1], B_prime.x])
C = B.orientnew('C', 'Axis', [q[2], B.z])

## --- Define Points and their velocities ---
pO = Point('O')
pO.set_vel(A, 0)

# R is the point in plane H that comes into contact with disk C.
pR = pO.locatenew('R', q[3]*A.x + q[4]*A.y)
pR.set_vel(A, pR.pos_from(pO).diff(t, A))
pR.set_vel(B, 0)

# C^ is the point in disk C that comes into contact with plane H.
pC_hat = pR.locatenew('C^', 0)
pC_hat.set_vel(C, 0)

# C* is the point at the center of disk C.
pCs = pC_hat.locatenew('C*', R*B.y)
pCs.set_vel(C, 0)
pCs.set_vel(B, 0)

# calculate velocities in A
pCs.v2pt_theory(pR, A, B)
pC_hat.v2pt_theory(pCs, A, C)

print("velocities of points R, C^, C* in rf A:")
print("v_R_A = {0}\nv_C^_A = {1}\nv_C*_A = {2}".format(
    pR.vel(A), pC_hat.vel(A), pCs.vel(A)))

## --- Expressions for generalized speeds u1, u2, u3, u4, u5 ---
u_expr = map(lambda x: dot(C.ang_vel_in(A), x), B)
u_expr += qd[3:]
kde = [u_i - u_ex for u_i, u_ex in zip(u, u_expr)]
kde_map = solve(kde, qd)
print("using the following kinematic eqs:\n{0}".format(msprint(kde)))

## --- Define forces on each point in the system ---
R_C_hat = Px*A.x + Py*A.y + Pz*A.z
R_Cs = -m*g*A.z
forces = [(pC_hat, R_C_hat), (pCs, R_Cs)]

## --- Calculate generalized active forces ---
partials = partial_velocities([pC_hat, pCs], u, A, kde_map)
F, _ = generalized_active_forces(partials, forces)
print("Generalized active forces:")
for i, f in enumerate(F, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))

# Now impose the condition that disk C is rolling without slipping
u_indep = u[:3]
u_dep = u[3:]
vc = map(lambda x: dot(pC_hat.vel(A), x), [A.x, A.y])
vc_map = solve(subs(vc, kde_map), u_dep)

partials_tilde = partial_velocities([pC_hat, pCs], u_indep, A, kde_map, vc_map)
F_tilde, _ = generalized_active_forces(partials_tilde, forces)
print("Nonholonomic generalized active forces:")
for i, f in enumerate(F_tilde, 1):
    print("F{0} = {1}".format(i, msprint(simplify(f))))

#print("\nUsing KanesMethod...")
#KM = KanesMethod(A, q, u, kde)
#fr, _ = KM.kanes_equations(forces, [])
#print("Generalized active forces:\n{0}".format(simplify(fr)))
#
#KM_tilde = KanesMethod(A, q, u_indep, kde,
#                       u_dependent=u_dep,
#                       velocity_constraints=vc)
#fr_tilde, _ = KM_tilde.kanes_equations(forces, [])
#print("Nonholonomic generalized active forces\n{0}".format(simplify(fr_tilde)))
"""Note: KanesMethod does not support a system that is defined independently of
u's. Also the KanesMethod calculation of the generalized active forces and
generalized inertia forces have the sign flipped.
"""
