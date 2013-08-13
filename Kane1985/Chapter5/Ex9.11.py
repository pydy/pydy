#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.11 from Kane 1985."""

from __future__ import division
from sympy import Dummy
from sympy import collect, expand, sin, cos, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, potential_energy
from util import generalized_active_forces_V


q = dynamicsymbols('q')
qd = dynamicsymbols('q', level=1)
u = dynamicsymbols('u')
L, L_prime, m, g = symbols('L L\' m g', real=True, positive=True)

# reference frames
N = ReferenceFrame('N')
# N.x points of the the plane of Figure P9.11, N.z points upward.
A = N.orientnew('A', 'axis', [q, N.x])


# define points
pO = Point('O') # point O is where the pendulum attaches to the ceiling
pP = pO.locatenew('P', -L * A.z) # mass center of the pendulum
pP.set_vel(N, pP.pos_from(pO).dt(N))

# kinematic differential equations
kde_map = {qd: u}

# forces
k = 5*m*g/L
r = (L_prime + L*sin(q))*N.y + (L - L*cos(q))*N.z
forces = [(pP, -m*g*N.z), (pP, -k*(r.magnitude() - L_prime)*r.normalize())]

partials = partial_velocities(zip(*forces)[0], [u], N, kde_map)
Fr, _ = generalized_active_forces(partials, forces)

# use a dummy symbol since series() does not work with dynamicsymbols
print('part a')
_q = Dummy('q')
terms = Fr[0].subs(q, _q).series(_q, n=4).removeO().subs(_q, q)
print('Using a series approximation of order 4:')
print('F1 ≈ {0}'.format(msprint(collect(terms, m*g*L))))

V = potential_energy([terms], [q], [u], kde_map)
print('V = {0}'.format(msprint(V)))
print('Setting C = 0, α1 = 0')
V = V.subs(dict(zip(symbols('C α1'), [0, 0])))
print('V = {0}'.format(msprint(collect(V, m*g*L))))

V_expected = m*g*L*(0*q + 3*q**2 + 0*q**3 + -7*q**4/8)
assert expand(V - V_expected) == 0


print('\npart b')
Fr_expected = m*g*L*(-6*q + 0*q**2 + 7*q**3/2)

print('Fr using V')
Fr_V = generalized_active_forces_V(V, [q], [u], kde_map)
print('F1_V = {0}'.format(msprint(collect(Fr_V[0], m*g*L))))
assert expand(Fr_V[0] - Fr_expected) == 0

print('Fr not using V, as calculated in part a')
print('F1 = {0}'.format(msprint(collect(terms, m*g*L))))
assert expand(terms - Fr_expected) == 0
