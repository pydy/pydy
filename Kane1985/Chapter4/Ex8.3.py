#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.3 from Kane 1985.
"""

from __future__ import division
#import os
#import sys
#sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
#                             "../../../sympy/sympy/physics/mechanics"))

from sympy import diff, factor, pi, solve, simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import partial_velocity
from sympy.physics.mechanics import MechanicsStrPrinter
from sympy.physics.mechanics import KanesMethod
#from kane import KanesMethod


def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)

def subs(x, *args, **kwargs):
    if not hasattr(x, 'subs'):
        if hasattr(x, '__iter__'):
            return map(lambda x: subs(x, *args, **kwargs), x)
    return x.subs(*args, **kwargs)

def generalized_active_forces(partial_velocities, resultants):
    assert len(partial_velocities) == len(resultants)
    Flist = []
    degree = len(partial_velocities[0])
    for r in range(degree):
        f = sum(dot(v_Pi[r], R_i)
                for v_Pi, R_i in zip(partial_velocities, resultants))
        Flist.append(factor(simplify(f)))
    return Flist

g, m, Px, Py, Pz, R, t = symbols('g m Px Py Pz R t')
q = dynamicsymbols('q1:6')
qd = dynamicsymbols('q1:6', level=1)
u = dynamicsymbols('u1:6')

## --- Define ReferenceFrames ---
A = ReferenceFrame('A')
B_prime = A.orientnew('B_prime', 'Axis', [q[1 - 1], A.z])
B = B_prime.orientnew('B', 'Axis', [pi/2 - q[2 - 1], B_prime.x])
C = B.orientnew('C', 'Axis', [q[3 - 1], B.z])

## --- Define Points and their velocities ---
pO = Point('O')
pO.set_vel(A, 0)

# R is the point in plane H that comes into contact with disk C.
pR = pO.locatenew('R', q[4 - 1]*A.x + q[5 - 1]*A.y)
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
u_expr += qd[4 - 1:]
kinematic_eqs = [u_i - u_ex for u_i, u_ex in zip(u, u_expr)]
soln = solve(kinematic_eqs, qd)
print("using the following kinematic eqs:\n{0}".format(msprint(kinematic_eqs)))

## --- Define forces on each point in the system ---
R_C_hat = Px*A.x + Py*A.y + Pz*A.z
R_Cs = -m*g*A.z
resultants = [R_C_hat, R_Cs]

## --- Calculate generalized active forces ---
vlist = subs([pC_hat.vel(A), pCs.vel(A)], soln)
v_r_Pi = partial_velocity(vlist, u, A)
F = generalized_active_forces(v_r_Pi, resultants)
print("Generalized active forces:")
for i, f in enumerate(F, 1):
    print("F{0} = {1}".format(i, msprint(f)))

# Now impose the condition that disk C is rolling without slipping
vel_constraint = subs(map(lambda x: dot(pC_hat.vel(A),x), [A.x, A.y]),
                           soln)
vel_cnstrnt_soln = solve(vel_constraint, u[4 - 1:])
vlist_tilde = subs(vlist, vel_cnstrnt_soln)
v_r_Pi_tilde = partial_velocity(vlist_tilde, u[:3], A)
F_tilde = generalized_active_forces(v_r_Pi_tilde, resultants)
print("Nonholonomic generalized active forces:")
for i, f in enumerate(F_tilde, 1):
    print("F{0} = {1}".format(i, msprint(f)))

print("\nUsing KanesMethod...")
KM = KanesMethod(A, q, u, kinematic_eqs)
FL = zip([pC_hat, pCs], resultants)
fr, _ = KM.kanes_equations(FL, [])
print("Generalized active forces:\n{0}".format(simplify(fr)))

u_indep = u[:3]
u_dep = u[4 - 1:]
KM_tilde = KanesMethod(A, q, u_indep, kinematic_eqs,
                       u_dependent=u_dep,
                       velocity_constraints=vel_constraint)
fr_tilde, _ = KM_tilde.kanes_equations(FL, [])
print("Nonholonomic generalized active forces\n{0}".format(simplify(fr_tilde)))
"""Note: KanesMethod does not support system that are defined independently of
u's. Also the KanesMethod calculation of the generalized active forces and
generalized inertia forces have the sign flipped.
"""
