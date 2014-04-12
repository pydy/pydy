#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 6.14 from Kane 1985."""

from __future__ import division
from sympy import Matrix, S
from sympy import integrate, pi, simplify, solve, symbols
from sympy.physics.mechanics import ReferenceFrame, Point
from sympy.physics.mechanics import cross, dot
from sympy.physics.mechanics import inertia, inertia_of_point_mass


def inertia_matrix(dyadic, rf):
    """Return the inertia matrix of a given dyadic for a specified
    reference frame.
    """
    return Matrix([[dot(dot(dyadic, i), j) for j in rf] for i in rf])


def integrate_v(integrand, rf, bounds):
    """Return the integral for a Vector integrand."""
    return sum(simplify(integrate(dot(integrand, n), bounds)) * n for n in rf)


def index_min(values):
    return min(range(len(values)), key=values.__getitem__)


m, a, b, c = symbols('m a b c', real=True, nonnegative=True)
x, y, r = symbols('x y r', real=True)
k, n = symbols('k n', real=True, positive=True)
N = ReferenceFrame('N')

# calculate right triangle density
V = b * c / 2
rho = m / V

# Kane 1985 lists scalars of inertia I_11, I_12, I_22 for a right triangular
# plate, but not I_3x.
pO = Point('O')
pCs = pO.locatenew('C*', b/3*N.x + c/3*N.y)
pP = pO.locatenew('P', x*N.x + y*N.y)
p = pP.pos_from(pCs)

I_3 = rho * integrate_v(integrate_v(cross(p, cross(N.z, p)),
                                  N, (x, 0, b*(1 - y/c))),
                      N, (y, 0, c))
print("I_3 = {0}".format(I_3))

# inertia for a right triangle given ReferenceFrame, height b, base c, mass
inertia_rt = lambda rf, b_, c_, m_: inertia(rf,
                m_*c_**2/18,
                m_*b_**2/18,
                dot(I_3, N.z),
                m_*b_*c_/36,
                dot(I_3, N.y),
                dot(I_3, N.x)).subs({m:m_, b:b_, c:c_})

theta = (30 + 90) * pi / 180
N2 = N.orientnew('N2', 'Axis', [theta, N.x])

# calculate the mass center of the assembly
# Point O is located at the right angle corner of triangle 1.
pCs_1 = pO.locatenew('C*_1', 1/S(3) * (a*N.y + b*N.x))
pCs_2 = pO.locatenew('C*_2', 1/S(3) * (a*N2.y + b*N2.x))
pBs = pO.locatenew('B*', 1/(2*m) * m * (pCs_1.pos_from(pO) +
                                        pCs_2.pos_from(pO)))
print("\nB* = {0}".format(pBs.pos_from(pO).express(N)))

# central inertia of each triangle
I1_C_Cs = inertia_rt(N, b, a, m)
I2_C_Cs = inertia_rt(N2, b, a, m)

# inertia of mass center of each triangle about Point B*
I1_Cs_Bs = inertia_of_point_mass(m, pCs_1.pos_from(pBs), N)
I2_Cs_Bs = inertia_of_point_mass(m, pCs_2.pos_from(pBs), N)

I_B_Bs = I1_C_Cs + I1_Cs_Bs + I2_C_Cs + I2_Cs_Bs
print("\nI_B_Bs = {0}".format(I_B_Bs.express(N)))

# central principal moments of inertia
evals = inertia_matrix(I_B_Bs, N).eigenvals()

print("\neigenvalues:")
for e in evals.iterkeys():
    print(e)

print("\nuse a/b = r")
evals_sub_a = [simplify(e.subs(a, r*b)) for e in evals.iterkeys()]
for e in evals_sub_a:
    print(e)

for r_val in [2, 1/S(2)]:
    print("\nfor r = {0}".format(r_val))
    evals_sub_r = [e.subs(r, r_val) for e in evals_sub_a]
    print("eigenvalues:")
    for e in evals_sub_r:
        print("{0} = {1}".format(e, e.n()))

    # substitute dummy values for m, b so that min will actually work
    min_index = index_min([e.subs({m : 1, b : 1}) for e in evals_sub_r])
    min_e = evals_sub_r[min_index]
    print("min: {0}".format(min_e))

    k_val = solve(min_e - 2*m*k**2, k)
    assert(len(k_val) == 1)
    print("k = {0}".format(k_val[0]))
    n_val = solve(k_val[0] - n*b, n)
    assert(len(n_val) == 1)
    print("n = {0}".format(n_val[0]))

print("\nResults in text: n = 1/3, sqrt(35 - sqrt(241))/24")
