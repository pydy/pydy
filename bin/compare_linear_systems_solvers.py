"""The purpose of this script is to compare the computational speed of the
various options for solving the dynamical differential equations linear system
for the accelerations, i.e M*u' = F -> u' = M^-1*F. Except for the
linear_sys_solver='sympy' option, PyDy executes a Python function for solving
Ax = b which has some overhead relative to doing the linear system solve in
C."""

import timeit

import numpy as np
import sympy as sm
import scipy as sp

from pydy import models
from pydy.codegen.ode_function_generators import (
    generate_ode_function, CythonODEFunctionGenerator)

sys = models.n_link_pendulum_on_cart(4, True, True)
constants = list(sm.ordered(sys.constants_symbols))
specifieds = list(sm.ordered(sys.specifieds_symbols))

t_val = 0.0
x_array = np.random.random(len(sys.states))
r_array = np.random.random(len(specifieds))
p_array = np.random.random(len(constants))

itr = 10000

rhs_sym = generate_ode_function(
    sys.eom_method.forcing,
    sys.coordinates,
    sys.speeds,
    constants,
    specifieds=specifieds,
    mass_matrix=sys.eom_method.mass_matrix,
    coordinate_derivatives=sm.Matrix(sys.speeds),
    linear_sys_solver='sympy',
    constants_arg_type='array',
    specifieds_arg_type='array',
    generator='cython')

time = timeit.timeit(lambda: rhs_sym(x_array, 0.0, r_array, p_array),
                     number=itr)
print('SymPy LU decomposition symbolic solve time: ', time)
sympy_solve_res = rhs_sym(x_array, 0.0, r_array, p_array)

g = CythonODEFunctionGenerator(
    sys.eom_method.forcing,
    sys.coordinates,
    sys.speeds,
    constants,
    specifieds=specifieds,
    mass_matrix=sys.eom_method.mass_matrix,
    coordinate_derivatives=sm.Matrix(sys.speeds),
    linear_sys_solver='numpy',
    constants_arg_type='array',
    specifieds_arg_type='array',
)
rhs_num = g.generate()

time = timeit.timeit(lambda: rhs_num(x_array, 0.0, r_array, p_array),
                     number=itr)
print('numpy.linalg.solve time: ', time)
numpy_linalg_solve_res = rhs_num(x_array, 0.0, r_array, p_array).copy()


def numpy_umath_linalg_solve(A, b):
    return np.squeeze(np.linalg._umath_linalg.solve(A, np.atleast_2d(b).T))


g.linear_sys_solver = numpy_umath_linalg_solve
time = timeit.timeit(lambda: rhs_num(x_array, 0.0, r_array, p_array),
                     number=itr)
print('numpy.linalg._umath_linalg.solve time: ', time)
numpy_umath_linalg_solve_res = rhs_num(x_array, 0.0, r_array, p_array).copy()


def scipy_linalg_solve(A, b):
    return sp.linalg.solve(A, b, lower=True, assume_a='pos',
                           check_finite=False, overwrite_a=True,
                           overwrite_b=True)


g.linear_sys_solver = scipy_linalg_solve
time = timeit.timeit(lambda: rhs_num(x_array, 0.0, r_array, p_array),
                     number=itr)
print('scipy.linalg.solve(..., assume_a="pos") time: ', time)
scipy_linalg_solve_res = rhs_num(x_array, 0.0, r_array, p_array).copy()


def scipy_linalg_lapack_dposv(A, b):
    _, x, _ = sp.linalg.lapack.dposv(A, b, lower=True, overwrite_a=True,
                                     overwrite_b=True)
    return np.squeeze(x)


g.linear_sys_solver = scipy_linalg_lapack_dposv
time = timeit.timeit(lambda: rhs_num(x_array, 0.0, r_array, p_array),
                     number=itr)
print('scipy.linalg.lapack.dposv time: ', time)
scipy_linalg_lapack_dposv_res = rhs_num(x_array, 0.0, r_array, p_array).copy()

np.testing.assert_allclose(sympy_solve_res,
                           numpy_linalg_solve_res)
np.testing.assert_allclose(numpy_linalg_solve_res,
                           numpy_umath_linalg_solve_res)
np.testing.assert_allclose(numpy_umath_linalg_solve_res,
                           scipy_linalg_solve_res)
np.testing.assert_allclose(scipy_linalg_solve_res,
                           scipy_linalg_lapack_dposv_res)
