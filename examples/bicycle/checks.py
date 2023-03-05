import numpy as np
import sympy as sm
import sympy.physics.mechanics as me
from pydy.codegen.ode_function_generators import CythonODEFunctionGenerator

from utils import (
    compare_numerical_arrays,
    evalf_with_symengine,
    evaluate_with_and_without_cse,
    evaluate_with_autowrap,
    formulate_equations_motion,
)


#print('Compare cse')

#compare_cse(F, args=[q3, q4, q5, q6, q7, q8, u3, u4, u5, u6, u7, u8, T4, T6,
                     #T7, rf, rr, d1, d2, d3, l1, l2, l3, l4, g, mc, md, me, mf,
                     #ic11, ic22, ic33, ic31, id11, id22, ie11, ie22, ie33,
                     #ie31, if11, if22])
#
#compare_cse(M, args=[q3, q4, q5, q6, q7, q8, u3, u4, u5, u6, u7, u8, T4, T6,
                     #T7, rf, rr, d1, d2, d3, l1, l2, l3, l4, g, mc, md, me, mf,
                     #ic11, ic22, ic33, ic31, id11, id22, ie11, ie22, ie33,
                     #ie31, if11, if22])

def check_ode_function_generator_variations(mass_matrix, forcing_vector, kane,
                                            symbolics, constant_substitutions,
                                            slow=False):

    print('Generating a right hand side function.')

    rhs_of_kin_diffs = sm.Matrix([kane.kindiffdict()[k]
                                  for k in kane.q.diff(symbolics['time'])])

    # min mass matrix csed with numpy linear solve
    g1 = CythonODEFunctionGenerator(forcing_vector,
                                    kane.q[:],
                                    kane.u[:],
                                    list(constant_substitutions.keys()),
                                    mass_matrix=mass_matrix,
                                    coordinate_derivatives=rhs_of_kin_diffs,
                                    specifieds=symbolics['specified quantities'],
                                    constants_arg_type='array',
                                    specifieds_arg_type='array',
                                    cse=True)
    print('Generating rhs')
    rhs1 = g1.generate()

    # min mass matrix csed with sympy linear solve
    g2 = CythonODEFunctionGenerator(forcing_vector,
                                    kane.q[:],
                                    kane.u[:],
                                    list(constant_substitutions.keys()),
                                    mass_matrix=mass_matrix,
                                    coordinate_derivatives=rhs_of_kin_diffs,
                                    specifieds=symbolics['specified quantities'],
                                    constants_arg_type='array',
                                    specifieds_arg_type='array',
                                    linear_sys_solver='sympy',
                                    cse=True)
    print('Generating rhs')
    rhs2 = g2.generate()

    np.testing.assert_allclose(rhs2, rhs1)

    if slow:
        # symbolically solved full rhs csed
        # This solves for the right hand side symbolically and generates code
        # from that expression.
        # NOTE : This takes like 12+ hours to compile.
        g3 = CythonODEFunctionGenerator(kane.rhs(),
                                        kane.q[:],
                                        kane.u[:],
                                        list(constant_substitutions.keys()),
                                        specifieds=symbolics['specified quantities'],
                                        constants_arg_type='array',
                                        specifieds_arg_type='array')
        print('Generating rhs')
        rhs3 = g3.generate()
        np.testing.assert_allclose(rhs3, rhs1)


def check_xreplace_against_exact(M_exact, F_exact, mass_matrix, forcing_vector,
                                 float_subs):

    print('Evaluating numerically with xreplace')
    M_from_xreplace = sm.matrix2numpy(mass_matrix.xreplace(float_subs),
                                      dtype=float)
    F_from_xreplace = sm.matrix2numpy(forcing_vector.xreplace(float_subs),
                                      dtype=float)

    compare_numerical_arrays(M_exact, M_from_xreplace,
                             name='Mass matrix from xreplace')
    compare_numerical_arrays(F_exact, F_from_xreplace,
                             name='Forcing vector from xreplace')


def check_cse(mass_matrix, forcing_vector, float_subs):

    # this runs with only the mass matrix but uses like 10+ GB of memory to
    # compile both of these show descrepancies in the steer equation u7
    M_no_cse, M_with_cse = evaluate_with_and_without_cse(mass_matrix,
                                                         float_subs)
    F_no_cse, F_with_cse = evaluate_with_and_without_cse(forcing_vector,
                                                         float_subs)
    compare_numerical_arrays(
        M_no_cse, M_with_cse,
        name='Mass matrix from CythonMatrixGenerator with cse')
    compare_numerical_arrays(
        F_no_cse, F_with_cse,
        name='Forcing vector from CythonMatrixGenerator with cse')


def check_autowrap_against_exact(M_exact, F_exact, M_sym, F_sym, float_subs):

    print('Evaluating with autowrap C')
    M_autowrap_c = evaluate_with_autowrap(M_sym, float_subs, language="C")
    F_autowrap_c = evaluate_with_autowrap(F_sym, float_subs, language="C")
    compare_numerical_arrays(M_exact, M_autowrap_c,
                             name='Mass matrix from autowrap C')
    compare_numerical_arrays(F_exact, F_autowrap_c,
                             name='Forcing vector from autowrap C')

    print('Evaluating with autowrap Fortran')
    M_autowrap_fortran = evaluate_with_autowrap(M_sym, float_subs,
                                                language="Fortran")
    F_autowrap_fortran = evaluate_with_autowrap(F_sym, float_subs,
                                                language="Fortran")
    compare_numerical_arrays(M_exact, M_autowrap_fortran,
                             name='Mass matrix from autowrap Fortran')
    compare_numerical_arrays(F_exact, F_autowrap_fortran,
                             name='Forcing vector from autowrap Fortran')


def check_kanes_equations(symbolics):
    """Verifies that the algorithm in ``KanesMethod`` for forming Kane's
    equations is equivalent to an independently written algorithm by evaluating
    the resulting equations to a high precision."""

    """This takes a while to run. I ran it once and got:
    Independent mass matrix is not correct. Here is the relative error:
    /home/moorepants/src/pydy/examples/bicycle/utils.py:39: RuntimeWarning: invalid value encountered in true_divide
    moore_output_from_sub = {k: v for k, v in zip(speed_deriv_names,
    [[-2. -2. -2. -2. -2. -2.]
    [-2. -2. -2. -2. -2. -2.]
    [-2. -2. -2. -2. -2. -2.]
    [-0. -0. -0.  0.  0.  0.]
    [-0. nan  0.  0.  0.  0.]
    [-0. nan  0. nan -0. nan]]
    Independent forcing vector is not correct. Here is the relative error:
    [[-2.]
    [-2.]
    [-2.]
    [ 0.]
    [-0.]
    [ 0.]]

    This is relative error and there is a value of -2 for all independent speed
    equations.

    TODO : Look into this!
    """

    print("Generating equations with sympy.physics.mechanics.KanesMethod.")
    kane = me.KanesMethod(
        symbolics['newtonian reference frame'],
        symbolics['independent generalized coordinates'],
        symbolics['independent generalized speeds'],
        kd_eqs=symbolics['kinematical differential equations'],
        q_dependent=symbolics['dependent generalized coordinates'],
        configuration_constraints=symbolics['holonomic constraints'],
        u_dependent=symbolics['dependent generalized speeds'],
        velocity_constraints=symbolics['nonholonomic constraints']
    )
    kane.kanes_equations(symbolics['bodies'], loads=symbolics['loads'])

    mass_matrix = kane.mass_matrix
    forcing_vector = kane.forcing.xreplace(kane.kindiffdict())

    print("Generating equations with independent method.")
    M, F = formulate_equations_motion(
        symbolics['newtonian reference frame'],
        symbolics['bodies'],
        symbolics['independent generalized coordinates'] +
        symbolics['dependent generalized coordinates'],
        symbolics['speed definitions'],
        symbolics['independent generalized speeds'],
        symbolics['dependent generalized speeds'],
        symbolics['nonholonomic constraints'],
        dict(symbolics['loads'])
    )

    print('Numerically evaluating the mass matrices and forcing vectors.')
    all_symbols = (
        symbolics['constants'] +
        symbolics['generalized coordinates'] +
        symbolics['generalized speeds'] +
        symbolics['specified quantities']
    )
    float_subs = dict(zip(all_symbols, np.random.random(len(all_symbols))))

    print('Evaluating KanesMethod mass matrix')
    M_sympy = evalf_with_symengine(mass_matrix, float_subs)
    print('Evaluating KanesMethod forcing vector')
    F_sympy = evalf_with_symengine(forcing_vector, float_subs)

    print('Evaluating independent method mass matrix')
    M_indep = evalf_with_symengine(M, float_subs)
    print('Evaluating independent method forcing vector')
    F_indep = evalf_with_symengine(F, float_subs)

    print('Comparing the numerical mass matrices and forcing vectors.')
    compare_numerical_arrays(M_indep, M_sympy,
                             name='Independent mass matrix')
    compare_numerical_arrays(F_indep, F_sympy,
                             name='Independent forcing vector')
