#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convenient utility functions for exercises in Chapter 4 of Kane 1985."""

from __future__ import division
from sympy.physics.mechanics import ReferenceFrame, Point, Particle, RigidBody
from sympy.physics.mechanics import cross, dot, Vector
from sympy.physics.mechanics import MechanicsStrPrinter
from sympy.physics.mechanics import inertia_of_point_mass


def msprint(expr):
    pr = MechanicsStrPrinter()
    return pr.doprint(expr)


def subs(x, *args, **kwargs):
    if x == 0:
        return x
    if not hasattr(x, 'subs'):
        if hasattr(x, '__iter__'):
            return map(lambda x: subs(x, *args, **kwargs), x)
    return x.subs(*args, **kwargs).doit()


class PartialVelocity(dict):
    def __init__(self, frame, ulist, *args, **kwargs):
        self._set_frame(frame)
        self._set_ulist(ulist)
        dict.__init__(self, *args, **kwargs)

    def _set_frame(self, f):
        if not isinstance(f, ReferenceFrame):
            raise TypeError(
                    '{0} is not an instance of ReferenceFrame'.format(f))
        self._frame = f

    def _set_ulist(self, u):
        if not isinstance(u, list):
            raise TypeError(
                    '{0} is not an instance of list'.format(f))
        self._ulist = u

    @property
    def frame(self):
        return self._frame

    @property
    def ulist(self):
        return self._ulist


def partial_velocities(system, generalized_speeds, frame,
                       kde_map=None, constraint_map=None, express_frame=None):
    partials = PartialVelocity(frame, generalized_speeds)
    if express_frame is None:
        express_frame = frame

    for p in system:
        if p in partials:
            continue
        if isinstance(p, Point):
            v = p.vel(frame)
        elif isinstance(p, ReferenceFrame):
            v = p.ang_vel_in(frame)
        if kde_map is not None:
            v = subs(v, kde_map)
        if constraint_map is not None:
            v = subs(v, constraint_map)

        v_r_p = {}
        for u in generalized_speeds:
            v_r_p[u] = Vector([]) if v == 0 else v.diff(u, express_frame)
        partials[p] = v_r_p
    return partials


def generalized_active_forces(partials, forces, uaux=None):
    # use the same frame used in calculating partial velocities
    ulist = partials.ulist

    if uaux is not None:
        uaux_zero = dict(zip(uaux, [0] * len(uaux)))

    Fr = [0] * len(ulist)
    for pf in forces:
        p = pf[0]
        f = pf[1]
        for i, u in enumerate(ulist):
            if partials[p][u] != 0 and f != 0:
                r = dot(partials[p][u], f)
                # if more than 2 args, 3rd is an integral function, where the
                # input is the integrand
                if len(pf) > 2:
                    r = pf[2](r)

                # auxilliary speeds have no effect on original active forces
                if uaux is not None and u not in uaux:
                    r = subs(r, uaux_zero)
                Fr[i] += r
    return Fr, ulist


def _calculate_T_star(rb, frame, kde_map, constraint_map, uaux):
    # get central inertia
    # I_S/O = I_S/S* + I_S*/O
    I = rb.inertia[0] - inertia_of_point_mass(rb.mass,
            rb.masscenter.pos_from(rb.inertia[1]), rb.frame)

    alpha = rb.frame.ang_acc_in(frame)
    omega = rb.frame.ang_vel_in(frame)
    if uaux is not None:
        # auxilliary speeds do not change alpha, omega
        # use doit() to evaluate terms such as
        # Derivative(0, t) to 0.
        uaux_zero = dict(zip(uaux, [0] * len(uaux)))
        alpha = subs(alpha, uaux_zero)
        omega = subs(omega, uaux_zero)
    if kde_map is not None:
        alpha = subs(alpha, kde_map)
        omega = subs(omega, kde_map)
    if constraint_map is not None:
        alpha = subs(alpha, constraint_map)
        omega = subs(omega, constraint_map)

    return -dot(alpha, I) - dot(cross(omega, I), omega)


def generalized_inertia_forces(partials, bodies,
                               kde_map=None, constraint_map=None,
                               uaux=None):
    # use the same frame used in calculating partial velocities
    ulist = partials.ulist
    frame = partials.frame

    if uaux is not None:
        uaux_zero = dict(zip(uaux, [0] * len(uaux)))

    Fr_star = [0] * len(ulist)
    for b in bodies:
        if isinstance(b, RigidBody):
            p = b.masscenter
            m = b.mass
        elif isinstance(b, Particle):
            p = b.point
            m = b.mass
        else:
            raise TypeError('{0} is not a RigidBody or Particle'.format(b))

        # get acceleration of point
        a = p.acc(frame)
        if uaux is not None:
            # auxilliary speeds do not change a
            a = subs(a, uaux_zero)
        if kde_map is not None:
            a = subs(a, kde_map)
        if constraint_map is not None:
            a = subs(a, constraint_map)


        # get T* for RigidBodys
        if isinstance(b, RigidBody):
            T_star = _calculate_T_star(b, frame, kde_map, constraint_map, uaux)

        for i, u in enumerate(ulist):
            force_term = 0
            torque_term = 0

            # inertia force term
            force_term = dot(partials[p][u], -m*a)

            # add inertia torque term for RigidBodys
            if isinstance(b, RigidBody):
                torque_term = dot(partials[b.frame][u], T_star)

            # auxilliary speeds have no effect on original inertia forces
            if uaux is not None and u not in uaux:
                force_term = subs(force_term, uaux_zero)
                torque_term = subs(torque_term, uaux_zero)

            Fr_star[i] += force_term + torque_term

    return Fr_star, ulist
