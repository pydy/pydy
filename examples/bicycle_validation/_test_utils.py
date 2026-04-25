import sympy as sm
import sympy.physics.mechanics as me

from pydy.models import n_link_pendulum_on_cart

from utils import (compare_cse, decompose_linear_parts, inv_of_3_by_3,
                   solve_for_qdots)

t = me.dynamicsymbols._t


def test_compare_cse():

    sys = n_link_pendulum_on_cart(n=5)

    F = sys.eom_method.forcing

    compare_cse(F)


def test_decompose_linear_parts():
    a, b, c, d, e, f, g, h, i = sm.symbols('a, b, c, d, e, f, g, h, i')
    u = me.dynamicsymbols('u0:5')
    uI = u[:2]
    uD = u[2:]

    expr1 = a*uI[0] + b*uI[1] + c*uD[0] + d*uD[1] + e*uD[2] + c**2
    expr2 = f*uI[0] + g*uI[1] + h*uD[0] + i*uD[1] + a*uD[2] + b**2
    expr3 = c*uI[0] + e*uI[1] + f*uD[0] + g*uD[1] + h*uD[2] + i**2

    motion_constraints = [expr1, expr2, expr3]

    A1, A2, B = decompose_linear_parts(motion_constraints, uI, uD)

    A_GuI_exp = sm.Matrix([[a, b], [f, g], [c, e]])
    A_GuD_exp = sm.Matrix([[c, d, e], [h, i, a], [f, g, h]])
    B_G_exp = sm.Matrix([c**2, b**2, i**2])

    assert A1 == A_GuI_exp
    assert A2 == A_GuD_exp
    assert B == B_G_exp

    A, B = decompose_linear_parts(motion_constraints, u)

    A_exp = sm.Matrix([[a, b, c, d, e],
                       [f, g, h, i, a],
                       [c, e, f, g, h]])

    assert A == A_exp
    assert B == B_G_exp

    A, B = decompose_linear_parts(motion_constraints, uD)
    exprs_for_uD_1 = A.LUsolve(-B)
    A1, A2, B2 = decompose_linear_parts(motion_constraints, uI, uD)
    exprs_for_uD_2 = A2.LUsolve(-A1*sm.Matrix(uI) - B2)

    assert sm.simplify(exprs_for_uD_1 - exprs_for_uD_2) == sm.zeros(3, 1)

    expr1 = uI[0] + uI[1] + a + b + c
    expr2 = uI[0] + e + f + g

    A, B = decompose_linear_parts([expr1, expr2], uI[0:2])
    res = A.LUsolve(-B)
    assert res[0] == -e - f - g
    assert res[1] == -a - b - c + e + f + g


def test_formulate_equations_of_motion():
    pass


def test_formulate_mass_matrix_form():
    pass


def test_generalized_active_forces():
    pass


def test_generalized_inertia_forces():
    pass


def test_inv_of_3_by_3():

    a, b, c, d, e, f, g, h, i = sm.symbols('a, b, c, d, e, f, g, h, i')

    A = sm.Matrix([[a, b, c], [d, e, f], [g, h, i]])

    assert sm.simplify(A.inv() - inv_of_3_by_3(A)) == sm.zeros(3, 3)


def test_solve_for_qdots():

    a = sm.symbols('a')
    q = me.dynamicsymbols('q0:3')
    u = me.dynamicsymbols('u0:3')

    u_def = {u[0]: q[0].diff(t),
             u[1]: q[1].diff(t)*sm.cos(q[2]) + a**2,
             u[2]: q[2].diff(t)*a*q[0]*q[1] + sm.sin(a)}

    res = solve_for_qdots(q, u_def)

    assert res[q[0].diff(t)] == u[0]
    assert sm.simplify(res[q[1].diff(t)] - (u[1] - a**2)/sm.cos(q[2])) == 0
    assert sm.simplify(res[q[2].diff(t)] - (u[2] - sm.sin(a))/a/q[0]/q[1]) == 0
