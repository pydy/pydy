#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 6.13 from Kane 1985
"""

from __future__ import division
from sympy import S, Matrix
from sympy import pi, acos
from sympy import simplify, sqrt, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import inertia, inertia_of_point_mass
from sympy.physics.mechanics import cross, dot

import numpy as np

def inertia_matrix(dyadic, rf):
    """Return the inertia matrix of a given dyadic for a specified
    reference frame.
    """
    return Matrix([[dot(dot(dyadic, i), j) for j in rf] for i in rf])

def convert_eigenvectors_matrix_vector(eigenvectors, rf):
    """Return a list of Vectors converted from a list of Matrices.
    rf is the implicit ReferenceFrame for the Matrix representation of the
    eigenvectors.
    """
    return [sum(simplify(v[0][i]).evalf() * n for i, n in enumerate(N))
            for v in eigenvectors]

def angle_between_vectors(a, b):
    """Return the minimum angle between two vectors. The angle returned for
    vectors a and -a is 0.
    """
    angle = (acos(dot(a, b)/(a.magnitude() * b.magnitude())) * 180 / pi).evalf()
    return min(angle, 180 - angle)

m = symbols('m', real=True, nonnegative=True)
m_val = 1
N = ReferenceFrame('N')
pO = Point('O')
pP = pO.locatenew('P', -3 * N.y)
pQ = pO.locatenew('Q', -4 * N.z)
pR = pO.locatenew('R', 2 * N.x)
points = [pO, pP, pQ, pR]

# center of mass of assembly
pCs = pO.locatenew('C*', sum(p.pos_from(pO) for p in points) / S(len(points)))
print(pCs.pos_from(pO))

I_C_Cs = sum(inertia_of_point_mass(m, p.pos_from(pCs), N) for p in points)
print("I_C_Cs = {0}".format(I_C_Cs))

# calculate the principal moments of inertia and the principal axes
M = inertia_matrix(I_C_Cs, N)

# use numpy to find eigenvalues/eigenvectors since sympy failed
# note that the eigenvlaues/eigenvectors are the
# prinicpal moments of inertia/principal axes
eigvals, eigvecs_np = np.linalg.eigh(np.matrix(M.subs(m, m_val).n().tolist()))
eigvecs = [sum(eigvecs_np[i, j] * n for i, n in enumerate(N)) for j in range(3)]

# get the minimum moment of inertia and its associated principal axis
e, v = min(zip(eigvals, eigvecs))

# I = m * k**2, where I is the moment of inertia,
# m is the mass of the body, k is the radius of gyration
k = sqrt(e / (4 * m_val))
print("\nradius of gyration, k = {0} m".format(k))

# calculate the angle between the associated principal axis and the line OP
# line OP is parallel to N.y
theta = angle_between_vectors(N.y, v)
print("\nangle between associated principal axis and line OP = {0}°".format(theta))
