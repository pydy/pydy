#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercises 8.1, 8.15 from Kane 1985.
"""

from __future__ import division
from collections import OrderedDict
from sympy import diff, factor, solve, simplify, symbols
from sympy import sin
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dot
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.mechanics import partial_velocity
from sympy.physics.mechanics import MechanicsStrPrinter

def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)

def subs(x, *args, **kwargs):
    if not hasattr(x, 'subs'):
        if hasattr(x, '__iter__'):
            return map(lambda x: subs(x, *args, **kwargs), x)
    return x.subs(*args, **kwargs)

class PartialVelocity(dict):
    def __init__(self, frame, *args, **kwargs):
        self._set_frame(frame)
        dict.__init__(self, *args, **kwargs)

    def _set_frame(self, f):
        if not isinstance(f, ReferenceFrame):
            raise TypeError(
                    '{0} is not an instance of ReferenceFrame'.format(f))
        self._frame = f

    @property
    def frame(self):
        return self._frame

def partial_velocities(system, generalized_speeds, frame,
                       kde_map=None, constraint_map=None, express_frame=None):
    partials = PartialVelocity(frame)
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

        v_r_p = OrderedDict()
        for u in generalized_speeds:
            v_r_p[u] = 0 if v == 0 else v.diff(u, express_frame)
        partials[p] = v_r_p
    return partials

def generalized_active_forces(partials, point_forces):
    # use the same frame used in calculating partial velocities
    v = partials.values()[0] # dict of partial velocities of the first item
    ulist = v.keys() # list of generalized speeds in case user wants it

    Fr = [0] * len(ulist)
    for ft in point_forces:
        p = ft[0]
        f = ft[1]
        for i, u in enumerate(ulist):
            if partials[p][u] != 0 and f != 0:
                r = dot(partials[p][u], f)
                if len(ft) == 2:
                    Fr[i] += r
                # if more than 2 args, 3rd is an integral function, where the
                # input is the integrand
                if len(ft) > 2:
                    Fr[i] += ft[2](r)
    return Fr, ulist

def generalized_inertia_forces(partials, point_masses, kde_map=None):
    # use the same frame used in calculating partial velocities
    v = partials.values()[0] # dict of partial velocities of the first item
    ulist = v.keys() # list of generalized speeds in case user wants it
    frame = partials.frame

    Fr_star = [0] * len(ulist)
    for p, m in point_masses:
        for i, u in enumerate(ulist):
            if partials[p][u] != 0 and m != 0:
                a = p.acc(frame)
                if kde_map is not None:
                    a = a.subs(kde_map)
                if a != 0:
                    Fr_star[i] += dot(partials[p][u], -m*a)
    return Fr_star, ulist

g, L, m1, m2, omega, t = symbols('g L m1 m2 omega t')
C, f1, f2 = symbols('C f1 f2')
q1, q2, q3 = dynamicsymbols('q1:4')
q1d, q2d, q3d = dynamicsymbols('q1:4', level=1)
u1, u2, u3 = dynamicsymbols('u1:4')

A = ReferenceFrame('A')
B = A.orientnew('B', 'Axis', [omega * t, A.y])
E = B.orientnew('E', 'Axis', [q3, B.z])

pO = Point('O')
pO.set_vel(A, 0)
pO.set_vel(B, 0)
pP1 = pO.locatenew('P1', q1 * B.x + q2 * B.y)
pP2 = pP1.locatenew('P2', L * E.x)
pP1.set_vel(E, 0)
pP1.set_vel(B, pP1.pos_from(pO).diff(t, B))
pP1.v1pt_theory(pO, A, B)
pP2.set_vel(E, 0)
pP2.v2pt_theory(pP1, A, E)

print("velocities of points P1, P2 in rf A:\nv_P1_A = {0}\nv_P2_A = {1}".format(
    pP1.vel(A), pP2.vel(A)))

# three sets of generalized speeds
u_s1 = [dot(pP1.vel(A), A.x), dot(pP1.vel(A), A.y), q3d]
u_s2 = [dot(pP1.vel(A), E.x), dot(pP1.vel(A), E.y), q3d]
u_s3 = [q1d, q2d, q3d]

# f1, f2 are forces the panes of glass exert on P1, P2 respectively
R1 = f1*B.z + C*E.x - m1*g*B.y
R2 = f2*B.z - C*E.x - m2*g*B.y

forces = [(pP1, R1), (pP2, R2)]
point_masses = [(pP1, m1), (pP2, m2)]
torques = []

ulist = [u1, u2, u3]
for uset in [u_s1, u_s2, u_s3]:
    print("\nFor generalized speeds:\n[u1, u2, u3] = {0}".format(msprint(uset)))
    # solve for u1, u2, u3 in terms of q1d, q2d, q3d and substitute
    kde = [u_i - u_expr for u_i, u_expr in zip(ulist, uset)]
    kde_map = solve(kde, [q1d, q2d, q3d])

    # include second derivatives in kde map
    for k, v in kde_map.items():
        kde_map[k.diff(t)] = v.diff(t)

    partials = partial_velocities([pP1, pP2], ulist, A, kde_map)
    Fr, _ = generalized_active_forces(partials, forces + torques)
    Fr_star, _ = generalized_inertia_forces(partials, point_masses, kde_map)
    print("Generalized active forces:")
    for i, f in enumerate(Fr, 1):
        print("F{0} = {1}".format(i, msprint(simplify(f))))
    print("Generalized inertia forces:")
    for i, f in enumerate(Fr_star, 1):
        sub_map = {}
        if uset == u_s1: # make the results easier to read
            if i == 1 or i == 3:
                sub_map = solve([u1 - u_s1[0]], [omega*q1*sin(omega*t)])
        print("F{0}* = {1}".format(i, msprint(simplify(f.subs(sub_map)))))

