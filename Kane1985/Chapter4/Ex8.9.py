#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.9 from Kane 1985.
"""

from __future__ import division
from collections import OrderedDict
from sympy import diff, solve, trigsimp, symbols, S
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import MechanicsStrPrinter


def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)

def subs(x, *args, **kwargs):
    if not hasattr(x, 'subs'):
        if hasattr(x, '__iter__'):
            return map(lambda x: subs(x, *args, **kwargs), x)
    return x.subs(*args, **kwargs)

def partial_velocities(system, generalized_speeds, frame,
                       kde_map=None, constraint_map=None, express_frame=None):
    partials = {}
    if express_frame is None:
        express_frame = frame

    for p in system:
        if isinstance(p, Point):
            v = p.vel(frame)
        elif isinstance(p, ReferenceFrame):
            v = p.ang_vel_in(frame)
        if kde_map is not None:
            v = v.subs(kde_map)
        if constraint_map is not None:
            v = v.subs(constraint_map)
        v_r_p = OrderedDict((u, v.diff(u, express_frame))
                            for u in generalized_speeds)
        partials[p] = v_r_p
    return partials

def generalized_active_forces(partials, force_pairs):
    # use the same frame used in calculating partial velocities
    v = partials.values()[0] # dict of partial velocities of the first item
    ulist = v.keys() # list of generalized speeds in case user wants it

    Flist = [0] * len(ulist)
    for p, f in force_pairs:
        for i, u in enumerate(ulist):
            if partials[p][u] and f:
                Flist[i] += dot(partials[p][u], f)
    return Flist, ulist

## --- Declare symbols ---
q1 = dynamicsymbols('q1')
q1d = dynamicsymbols('q1', level=1)
u1 = dynamicsymbols('u1')
g, m, r, theta = symbols('g m r theta')

print("part a")
# --- ReferenceFrames ---
N = ReferenceFrame('N')
C = N.orientnew('C', 'Axis', [q1, N.x])
R = C # R is fixed relative to C
A = C.orientnew('A', 'Axis', [-theta, R.x])
B = C.orientnew('B', 'Axis', [theta, R.x])

# --- Define points and their velocities ---
pO = Point('O') # Point O is at the center of the circular wire
pA = pO.locatenew('A', -r * A.z)
pB = pO.locatenew('B', -r * B.z)
pR_star = pO.locatenew('R*', 1/S(2) * (pA.pos_from(pO) + pB.pos_from(pO)))

pO.set_vel(N, 0)
pO.set_vel(C, 0)
for p in [pA, pB, pR_star]:
    p.set_vel(C, 0)
    p.v1pt_theory(pO, N, C)

## --- Define kinematic differential equations ---
kde = [u1 - q1d]
kde_map = solve(kde, [q1d])

## --- Define contact/distance forces ---
forces = [(pR_star, -m*g*N.z)]
torques = []

def print_fr(forces, ulist):
    print("Generalized active forces:")
    partials = partial_velocities(zip(*forces + torques)[0], ulist, N, kde_map)
    Fr, _ = generalized_active_forces(partials, forces + torques)
    for i, f in enumerate(Fr, 1):
        print("F{0} = {1}".format(i, msprint(trigsimp(f)))) 
print_fr(forces, [u1])

print("\npart b")
# define new symbols and incorporate new friction forces in fr
s, TA, TB = symbols('s TA TB')
forces += [(pA, -s*TA*A.y), (pB, -s*TB*B.y)]

print_fr(forces, [u1])

print("\npart b")
# define u2, u3, and radial forces and incorporate in fr
u2, u3 = dynamicsymbols('u2 u3')
NA, NB = symbols('NA NB')

# must redefine velocities of A, B, R* since R is no longer fixed in C
for p in [pA, pB, pR_star]:
    p.set_vel(R, 0)
pR_star.set_vel(C, u2*C.y + u3*C.z)
for p in [pA, pB]:
    p.v1pt_theory(pR_star, C, R)
for p in [pA, pB, pR_star]:
    p.v1pt_theory(pO, N, C)
forces += [(pA, NA*A.z), (pB, NB*B.z)]
print_fr(forces, [u1, u2, u3])

