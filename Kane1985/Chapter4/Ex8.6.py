#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 8.6 from Kane 1985.
"""

from __future__ import division
from collections import OrderedDict
from sympy import diff, solve, simplify, symbols
from sympy import pi, sin, cos
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
# Define the system with 6 generalized speeds as follows:
q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')
L1, L2, L3, L4 = symbols('L1:5')
g, m1, m2, t = symbols('g m1 m2 t')

# --- ReferenceFrames ---
A = ReferenceFrame('A')

# --- Define Points and set their velocities ---
pO = Point('O')
pO.set_vel(A, 0)
pP1 = pO.locatenew('P1', L1*(cos(q1)*A.x + sin(q1)*A.y))
pP1.set_vel(A, pP1.pos_from(pO).diff(t, A))
pP2 = pP1.locatenew('P2', L2*(cos(q2)*A.x + sin(q2)*A.y))
pP2.set_vel(A, pP2.pos_from(pO).diff(t, A))

## --- configuration constraints ---
cc = [L1*cos(q1) + L2*cos(q2) - L3*cos(q3),
      L1*sin(q1) + L2*sin(q2) - L3*sin(q3) - L4]

## --- Define kinematic differential equations/pseudo-generalized speeds ---
kde = [u1 - q1d, u2 - q2d, u3 - q3d]
kde_map = solve(kde, [q1d, q2d, q3d])

# --- velocity constraints ---
vc = [c.diff(t) for c in cc]
vc_map = solve(subs(vc, kde_map), [u2, u3])

## --- Define gravitational forces ---
forces = [(pP1, m1*g*A.x), (pP2, m2*g*A.x)]

partials = partial_velocities([pP1, pP2], [u1], A, kde_map, vc_map)
Fr, _ = generalized_active_forces(partials, forces)
print("Generalized active forces:")
for i, f in enumerate(Fr, 1):
    print("F{0}_tilde = {1}".format(i, simplify(f)))

