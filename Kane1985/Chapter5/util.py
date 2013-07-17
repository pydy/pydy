#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convenient utility functions for exercises in Chapter 5 of Kane 1985."""

from __future__ import division
from sympy import Dummy, Matrix
from sympy import diff, expand, expand_trig, integrate, solve, symbols
from sympy import sin, cos, tan, trigsimp
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
        else:
            raise TypeError('A Point or ReferenceFrame must be supplied.')
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
        p = pf[0] # first arg is point/rf
        f = pf[1] # second arg is force/torque
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


def _equivalent_derivatives(dV_dq_list, q):
    dV_eq = []
    for r in range(len(q)):
        for s in range(r + 1, len(q)):
            dV_eq.append(dV_dq_list[r].diff(q[s]) - dV_dq_list[s].diff(q[r]))
    return dV_eq


def _f_variables(Fr, q, dV_eq, dV_dq):
    Fr_qi_only = []
    non_arg = set()
    for i, fr in enumerate(Fr):
        dfrdqi = [j for j, x in enumerate(q) if fr.diff(x) != 0]
        # If generalized force is only a function of one generalized coordinate
        # save the indices of force, coordinate.
        if len(dfrdqi) == 1:
            Fr_qi_only.append((i, dfrdqi[0]))

    for fr_idx, qi_idx in Fr_qi_only:
        # If Fr = -∂V/∂qi, then fs-p is independent of qi.
        if Fr[fr_idx] - dV_eq[fr_idx] == dV_dq[qi_idx]:
            non_arg.add(q[qi_idx])
    return sorted(list(set(q) - non_arg)) + [symbols('t')]


def potential_energy(Fr, q, u, kde_map, vc_map=None):
    if vc_map is not None:
        u += sorted(vc_map.keys())

    m = len(u) - len(Fr)
    dV_dq = symbols('∂V/∂q1:{0}'.format(len(q) + 1))
    dV_eq = Matrix(Fr).T
    W_sr = Matrix([map(lambda x: v.diff(x), u)
                   for k, v in sorted(kde_map.iteritems())])
    if vc_map is not None:
        A_kr = Matrix([map(lambda x: diff(v, x), u[:len(Fr)])
                       for k, v in sorted(vc_map.iteritems())])
    else:
        A_kr = Matrix.zeros(m, len(Fr))

    for s in range(W_sr.shape[0]):
        dV_eq += dV_dq[s] * (W_sr[s, :len(Fr)] + W_sr[s, len(Fr):]*A_kr)

    if vc_map is not None:
        f = map(lambda x: x(*_f_variables(Fr, q, dV_eq, dV_dq)),
                symbols('f1:{0}'.format(m + 1)))
        dV_eq = subs(dV_eq, dict(zip(dV_dq[-m:], f)))
        dV_dq = dV_dq[:-m]

    dV_dq_map = solve(dV_eq, dV_dq)
    dV_dq_list = map(lambda x: dV_dq_map[x], dV_dq)

    if vc_map is None:
        print('Checking ∂/∂qr(∂V/∂qs) = ∂/∂qs(∂V/∂qr) for all r, s '
              '= 1, ..., n.')
        dV_eq = _equivalent_derivatives(dV_dq_list, q)
        if dV_eq != [0] * len(q):
            rs = [(r, s) for r in range(len(q)) for s in range(r + 1, len(q))]
            for (r, s), x in zip(rs, dV_eq):
                if x != 0:
                    print(('∂/∂q{0}(∂V/∂q{1}) != ∂/∂q{1}(∂V/∂q{0}). ' +
                           'V does NOT exist.').format(r + 1, s + 1))
                    print('∂/∂q{0}(∂V/∂q{1}) = {2}'.format(
                            r + 1, s + 1, dV_dq_list[r].diff(q[s])))
                    print('∂/∂q{1}(∂V/∂q{0}) = {2}'.format(
                            r + 1, s + 1, dV_dq_list[s].diff(q[r])))
                    break
            return None
    else:
        dV_dq_list += f
        # Unable to take diff of 'fm.diff(qs)', replace with dummy symbols.
        dfdq = [Dummy('∂f{0}/∂q{1}'.format(i + 1, j + 1))
                for i in range(len(f)) for j in range(len(q))]
        dfdq_replace = lambda x: reduce(
                lambda y, z: y.replace(z[0], z[1]),
                zip([fm.diff(qs) for fm in f for qs in q], dfdq),
                x)
        dV_eq = map(dfdq_replace,
                    _equivalent_derivatives(dV_dq_list, q))

        X = Matrix(dfdq)
        Z = Matrix([map(lambda x: diff(dV_eqi, x), dfdq)
                    for dV_eqi in dV_eq])
        if Z.rank() == len(q) * (len(q) - 1) / 2:
            print('ρ == n(n - 1)/2')
            print('V may exist but cannot be found by this procedure.')
            return None

        Y = expand(Z*X - Matrix(dV_eq))
        ZI_rref, _ = Matrix.hstack(Z, Matrix.eye(Z.shape[0])).rref()
        # E is the matrix of elementary row operations that gives rref(Z).
        E = ZI_rref[:, Z.shape[1]:]
        f_eq = (E * Y)[Z.rank():]
        f_map = solve(f_eq, f)
        dV_dq_list = map(trigsimp, (subs(dV_dq_list, f_map)))

    alpha = symbols('α1:{0}'.format(len(q) + 1))
    V, zeta = symbols('C ζ')
    q_alpha = zip(q, alpha)
    for i, dV_dqr in enumerate(dV_dq_list):
        integrand = dV_dqr.subs(dict(q_alpha[i + 1:])).subs(q[i], zeta)
        V += integrate(expand_trig(integrand), (zeta, alpha[i], q[i]))
    return V
