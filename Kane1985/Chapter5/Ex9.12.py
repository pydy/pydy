#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 9.12 from Kane 1985.
Answer does not match text.
"""

from __future__ import division
from sympy import Dummy
from sympy import expand, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import dynamicsymbols
from util import msprint, subs, partial_velocities
from util import generalized_active_forces, potential_energy


q1, q2, q3, q4, q5, q6 = q = dynamicsymbols('q1:7')
u1, u2, u3, u4, u5, u6 = u = dynamicsymbols('u1:7')
# L' is the natural length of the springs
a, k, L_prime = symbols('a k L\'', real=True, positive=True)

# reference frames
X = ReferenceFrame('X')
C = X.orientnew('C', 'body', [q4, q5, q6], 'xyz')

# define points
pO = Point('O') # point O is fixed in X
pC_star = pO.locatenew('C*', a*(q1*X.x + q2*X.y + q3*X.z))
# define points of the cube connected to springs
pC1 = pC_star.locatenew('C1', a*(C.x + C.y - C.z))
pC2 = pC_star.locatenew('C2', a*(C.y + C.z - C.x))
pC3 = pC_star.locatenew('C3', a*(C.z + C.x - C.y))
# define fixed spring points
pk1 = pO.locatenew('k1', L_prime * X.x + a*(X.x + X.y - X.z))
pk2 = pO.locatenew('k2', L_prime * X.y + a*(X.y + X.z - X.x))
pk3 = pO.locatenew('k3', L_prime * X.z + a*(X.z + X.x - X.y))

pC_star.set_vel(X, pC_star.pos_from(pO).dt(X))
pC1.v2pt_theory(pC_star, X, C)
pC2.v2pt_theory(pC_star, X, C)
pC3.v2pt_theory(pC_star, X, C)

# kinematic differential equations
kde_map = dict(zip(map(lambda x: x.diff(), q), u))

# forces
x1 = pC1.pos_from(pk1)
x2 = pC2.pos_from(pk2)
x3 = pC3.pos_from(pk3)
forces = [(pC1, -k*(x1.magnitude() - L_prime)*x1.normalize()),
          (pC2, -k*(x2.magnitude() - L_prime)*x2.normalize()),
          (pC3, -k*(x3.magnitude() - L_prime)*x3.normalize())]

partials = partial_velocities(zip(*forces)[0], u, X, kde_map)
Fr, _ = generalized_active_forces(partials, forces)
print('generalized active forces')
for i, fr in enumerate(Fr, 1):
    print('\nF{0} = {1}'.format(i, msprint(fr)))

# use a dummy symbol since series() does not work with dynamicsymbols
_q = Dummy('q')
series_exp = (lambda x, qi, n_:
                    x.subs(qi, _q).series(_q, n=n_).removeO().subs(_q, qi))

# remove all terms order 3 or higher in qi
Fr_series = [reduce(lambda x, y: series_exp(x, y, 3), q, fr)
             for fr in Fr]
print('\nseries expansion of generalized active forces')
for i, fr in enumerate(Fr_series, 1):
    print('\nF{0} = {1}'.format(i, msprint(fr)))

V = potential_energy(Fr_series, q, u, kde_map)
print('\nV = {0}'.format(msprint(V)))
print('Setting C = 0, α1, α2, α3, α4, α5, α6 = 0')
V = V.subs(dict(zip(symbols('C α1 α2 α3 α4 α5 α6'), [0] * 7)))
print('V = {0}'.format(msprint(V)))

V_expected = k*a**2/2*((q1 - q5 - q6)**2 + (q2 - q6 - q4)**2 +
                       (q3 - q4 - q5)**2)
assert expand(V - V_expected) == 0
