#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.10 from Kane 1985.
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
            if partials[p][u] != 0 and f != 0:
                Flist[i] += dot(partials[p][u], f)
    return Flist, ulist

## --- Declare symbols ---
u1, u2, u3, u4, u5 = dynamicsymbols('u1:6')
TB, TE = symbols('TB TE')

# --- ReferenceFrames ---
A = ReferenceFrame('A')
B = ReferenceFrame('B')
CD = ReferenceFrame('CD')
CD_prime = ReferenceFrame('CD_prime')
E = ReferenceFrame('E')
F = ReferenceFrame('F')

# --- Define angular velocities of reference frames ---
B.set_ang_vel(A, u1 * A.x)
B.set_ang_vel(F, u2 * A.x)
CD.set_ang_vel(F, u3 * F.y)
CD_prime.set_ang_vel(F, u4 * -F.y)
E.set_ang_vel(F, u5 * -A.x)

# --- define velocity constraints ---
teeth = dict([('A', 60), ('B', 30), ('C', 30),
              ('D', 61), ('E', 20)])
vc = [u2*teeth['B'] - u3*teeth['C'],  # w_B_F * r_B = w_CD_F * r_C
      u2*teeth['B'] - u4*teeth['C'],  # w_B_F * r_B = w_CD'_F * r_C
      u5*teeth['E'] - u3*teeth['D'],  # w_E_F * r_E = w_CD_F * r_D
      u5*teeth['E'] - u4*teeth['D'],  # w_E_F * r_E = w_CD'_F * r_D;
      (-u1 + u2)*teeth['A'] - u3*teeth['D'],  # w_A_F * r_A = w_CD_F * r_D;
      (-u1 + u2)*teeth['A'] - u4*teeth['D']]  # w_A_F * r_A = w_CD'_F * r_D;
vc_map = solve(vc, [u2, u3, u4, u5])

## --- Define torques ---
forces = []
torques = [(B, TB*A.x), (E, TE*A.x)]

partials = partial_velocities([B, E], [u1], A, constraint_map = vc_map)
Fr, _ = generalized_active_forces(partials, forces + torques)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0} = {1}".format(i, msprint(trigsimp(f)))) 
