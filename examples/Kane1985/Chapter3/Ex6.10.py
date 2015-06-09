#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 6.10 from Kane 1985."""

from __future__ import division
from sympy import Matrix
from sympy import pi, acos
from sympy import simplify, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia, inertia_of_point_mass
from sympy.physics.mechanics import dot


def inertia_matrix(dyadic, rf):
    """Return the inertia matrix of a given dyadic for a specified
    reference frame.
    """
    return Matrix([[dot(dot(dyadic, i), j) for j in rf] for i in rf])


def angle_between_vectors(a, b):
    """Return the minimum angle between two vectors. The angle returned for
    vectors a and -a is 0.
    """
    angle = (acos(dot(a, b)/(a.magnitude() * b.magnitude())) *
             180 / pi).evalf()
    return min(angle, 180 - angle)


m, m_R, m_C, rho, r = symbols('m m_R m_C rho r', real=True, nonnegative=True)

N = ReferenceFrame('N')
pA = Point('A')
pPs = pA.locatenew('P*', 3*r*N.x - 2*r*N.y)

m_R = rho * 24 * r**2
m_C = rho * pi * r**2
m = m_R - m_C

I_Cs_A = inertia_of_point_mass(m, pPs.pos_from(pA), N)
I_C_Cs = inertia(N, m_R*(4*r)**2/12 - m_C*r**2/4,
                    m_R*(6*r)**2/12 - m_C*r**2/4,
                    m_R*((4*r)**2+(6*r)**2)/12 - m_C*r**2/2)

I_C_A = I_C_Cs + I_Cs_A
print("\nI_C_rel_A = {0}".format(I_C_A))

# Eigenvectors of I_C_A are the parallel to the principal axis for point A
# of Body C.
evecs_m = [triple[2]
           for triple in inertia_matrix(I_C_A, N).eigenvects()]

# Convert eigenvectors from Matrix type to Vector type.
evecs = [sum(simplify(v[0][i]).evalf() * n for i, n in enumerate(N))
         for v in evecs_m]

# N.x is parallel to line AB
print("\nVectors parallel to the principal axis for point A of Body C and the" +
      "\ncorresponding angle between the principal axis and line AB (degrees):")
for v in evecs:
    print("{0}\t{1}".format(v, angle_between_vectors(N.x, v)))
