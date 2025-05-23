#!/usr/bin/env python

from collections import OrderedDict

import numpy as np
import sympy as sm
import sympy.physics.mechanics as me
from sympy.utilities.autowrap import autowrap
from pydy.codegen.cython_code import CythonMatrixGenerator
from dtk import bicycle
import symengine as se

TIME = me.dynamicsymbols._t


class ReferenceFrame(me.ReferenceFrame):
    """Subclass that enforces the desired unit vector indice style."""

    def __init__(self, *args, **kwargs):

        kwargs.pop('indices', None)
        kwargs.pop('latexs', None)

        lab = args[0].lower()
        tex = r'\hat{{{}}}_{}'

        super(ReferenceFrame, self).__init__(*args, indices=('1', '2', '3'),
                                             latexs=(tex.format(lab, '1'),
                                                     tex.format(lab, '2'),
                                                     tex.format(lab, '3')),
                                             **kwargs)


def create_basu_output_from_moore_output(symbolics, kane, g, xd_from_sub,
                                         xd_from_gen, substitutions,
                                         moore_input, bp):

    print("Generating output dictionary.")
    # convert the outputs from my model to the Basu-Mandal coordinates
    speed_deriv_names = [str(speed)[:-3] + 'p' for speed in kane.u[:]]

    moore_output_from_sub = {k: v for k, v in zip(speed_deriv_names,
                                                list(xd_from_sub))}
    moore_output_from_gen = {k: v for k, v in zip(speed_deriv_names,
                                                list(xd_from_gen)[g.num_coordinates:])}
    u5 = symbolics['dependent generalized speeds'][1]
    u6 = symbolics['independent generalized speeds'][1]
    t = symbolics['time']
    acc_subs_sub = {u5.diff(t): moore_output_from_sub['u5p'],
                    u6.diff(t): moore_output_from_sub['u6p']}

    acc_subs_gen = {u5.diff(t): moore_output_from_gen['u5p'],
                    u6.diff(t): moore_output_from_gen['u6p']}

    u1p_sym = sm.Function('u1')(t).diff(t)
    u2p_sym = sm.Function('u2')(t).diff(t)
    u1p = symbolics['extra definitions'][u1p_sym]
    u2p = symbolics['extra definitions'][u2p_sym]
    moore_output_from_sub['u1p'] = u1p.xreplace(acc_subs_sub).xreplace(substitutions)
    moore_output_from_sub['u2p'] = u2p.xreplace(acc_subs_sub).xreplace(substitutions)

    moore_output_from_gen['u1p'] = u1p.xreplace(acc_subs_gen).xreplace(substitutions)
    moore_output_from_gen['u2p'] = u2p.xreplace(acc_subs_gen).xreplace(substitutions)

    moore_output_from_sub.update(moore_input)
    moore_output_from_gen.update(moore_input)

    moore_output_basu_from_sub = bicycle.moore_to_basu(moore_output_from_sub,
                                                       bp['rR'], bp['lam'])
    moore_output_basu_from_gen = bicycle.moore_to_basu(moore_output_from_gen,
                                                       bp['rR'], bp['lam'])

    return moore_output_basu_from_sub, moore_output_basu_from_gen


def create_moore_input_from_basu_input(symbolics):

    print('Loading numerical input parameters.')

    # these are the Benchmark bicycle [Meijaard 2007] parameters.
    bp = bicycle.benchmark_parameters()
    # these are the Benchmark bicycle [Meijaard 2007] parameters reexpressed in
    # the [Moore 2012] definition.
    mp = bicycle.benchmark_to_moore(bp)

    # load the input values specified in Table 1 of [Basu-Mandal 2007]
    basu_input = bicycle.basu_table_one_input()

    # convert the Basu-Mandal values to my coordinates and speeds
    moore_input = bicycle.basu_to_moore_input(basu_input, bp['rR'], bp['lam'])

    constants_name_map = {sym.name: sym for sym in symbolics['constants']}
    time_varying_name_map = {s.name: s for s in
                            symbolics['generalized coordinates'] +
                            symbolics['generalized speeds'] +
                            symbolics['specified quantities']}

    # build dictionaries that map the Moore symbolic parameters to the
    # converted Basu-Mandal values
    constant_substitutions = OrderedDict()
    for k, v in mp.items():
        try:
            constant_substitutions[constants_name_map[k]] = v
        except KeyError:
            print('{} not added to sub dict.'.format(k))

    dynamic_substitutions = {}
    for k, v in moore_input.items():
        try:
            dynamic_substitutions[time_varying_name_map[k]] = v
        except KeyError:
            print('{} not added to sub dict.'.format(k))
        # TODO : try this to ensure we are using 0.0 instead of other tiny
        # floats.  As some of the converted Basu-Mandal numbers could have
        # floating point round off.
        #else:
            #if abs(dynamic_substitutions[k]) < 1e-14:
                #dynamic_substitutions[k] = 0.0

    # set all specifieds to zero
    specified_subs = {ri: 0.0 for ri in  symbolics['specified quantities']}

    constants_substituions, dynamic_substitution, specified_subs = \
        create_symbol_value_map(constants_name_map, time_varying_name_map)

    substitutions = specified_subs.copy()
    substitutions.update(constant_substitutions)
    substitutions.update(dynamic_substitutions)

    return (substitutions, constants_substituions, dynamic_substitution,
            specified_subs, moore_input, bp)



def compare_to_basu_values(exact_vals, float_vals):
    """Prints 16 digit comparisons and tests with
    numpy.testing.assert_allclose() to validate a match or not.

    Parameters
    ==========
    exact_vals : dictionary
        Maps Basu-Mandal variable strings to values computed from fixed
        precision methods.
    float_vals : dictionary
        Maps Basu-Mandal variable strings to values computed from floating
        point methods.

    """

    print('Comparing numerical results to the Basu-Mandal values:')

    fail_msg = 'Failed: {}:\n  Expected: {:1.16f}\n    Actual: {:1.16f}'
    matc_msg = 'Matched: {}:\n  Expected: {:1.16f}\n    Actual: {:1.16f}'

    # load values from Basu-Mandal 2007 Table 1
    basu_output = bicycle.basu_table_one_output()

    for result, typ in zip([exact_vals, float_vals],
                           ['Results from Symengine evaluation',
                            'Results from CythonODEGenerator']):
        print(typ)
        for k, bv in basu_output.items():
            try:
                mv = float(result[k])
            except KeyError:
                print('{} was not checked.'.format(k))
            else:
                try:
                    np.testing.assert_allclose(bv, mv)
                except AssertionError:
                    print(fail_msg.format(k, bv, mv))
                else:
                    print(matc_msg.format(k, bv, mv))


def compare_numerical_arrays(actual, expected, name='Actual'):
    try:
        np.testing.assert_allclose(actual, expected)
    except AssertionError:
        print('{} is not correct. Here is the relative error:'.format(name))
        print((actual - expected)/expected)
    else:
        print('{} is correct to machine precision'.format(name))


def evaluate_with_autowrap(expr, float_subs, language='C', tmp_dir=None):
    """Returns the numerical evaluations of a single symbolic matrix expression
    using ``autowrap``.

    Parameters
    ==========
    expr : sympy.Matrix, shape(n, m)
        A matrix of expressions.
    float_subs : dictionary
        Maps all the symbols in ``expr`` to floating point numbers.
    language : string
        "C" or "Fortran"
    tmp_dir : string
        Path to a directory to store the generated files.

    Returns
    =======
    res : ndarray, shape(n, m)

    """

    # autowrap can't handle the Function()(t), so subs in symbols.
    func_subs = {}
    new_float_subs = {}
    for k, v in float_subs.items():
        if '(t)' in str(k):
            s = sm.Symbol(k.name[-3:])
            func_subs[k] = s
            new_float_subs[s] = v
        else:
            new_float_subs[k] = v

    if language == 'C':
        backend = 'cython'
    elif language == 'Fortran':
        language = 'F95'
        backend = 'f2py'

    eval_expr = autowrap(expr.xreplace(func_subs),
                         language=language,
                         backend=backend,
                         tempdir=tmp_dir,
                         args=tuple(new_float_subs.keys()))
    res = eval_expr(*new_float_subs.values())

    return res


def evaluate_with_and_without_cse(expr, float_subs, tmp_dir=None):
    """Returns the numerical evaluations of a single symbolic matrix expression
    using ``CythonMatrixGenerator`` with and without using common subexpression
    elimination.

    Parameters
    ==========
    expr : sympy.Matrix, shape(n, m)
        A matrix of expressions.
    float_subs : dictionary
        Maps all the symbols in ``expr`` to floating point numbers.

    Returns
    =======
    res_no_cse : ndarray, shape(n, m)
    res_with_cse : ndarray, shape(n, m)

    """

    args = []
    vals = []
    for k, v in float_subs.items():
        args.append(k)
        vals.append(v)
    vals = np.array(vals)
    num_rows, num_cols = expr.shape

    gen_no_cse = CythonMatrixGenerator([args], [expr], cse=False)
    eval_no_cse = gen_no_cse.compile(tmp_dir=tmp_dir, verbose=True)
    res_no_cse = np.empty(num_rows*num_cols)
    eval_no_cse(vals, res_no_cse)

    gen_with_cse = CythonMatrixGenerator([args], [expr], cse=True)
    eval_with_cse = gen_with_cse.compile(tmp_dir=tmp_dir, verbose=True)
    res_with_cse = np.empty(num_rows*num_cols)
    eval_with_cse(vals, res_with_cse)

    return res_no_cse, res_with_cse


def evalf_with_symengine(sympy_expr, float_subs):
    """Returns the numerical evaluation of the expression using floats with X
    digits of precision using symengine for speed (SymPy's evalf() takes hours
    to compute).

    Parameters
    ==========
    expr : sympy.Matrix, shape(n, m)
        A matrix of expressions.
    float_subs : dictionary
        Maps all the symbols in ``expr`` to floating point numbers.

    Returns
    =======
    res : ndarray, shape(n, m)
        Array of double precision floats.

    """

    symengine_expr = se.sympify(sympy_expr)

    rational_subs = {}
    for k, v in float_subs.items():
        # 17 because it is needed to represent a double
        # https://stackoverflow.com/questions/6118231/why-do-i-need-17-significant-digits-and-not-16-to-represent-a-double/
        rational_subs[k] = se.Integer(int(v*10**17))/10**17
        #rational_subs[k] = se.Integer(int(v*2**52))/2**52

    symengine_expr_with_rationals = symengine_expr.xreplace(rational_subs)
    # n(number_of_bits, real=True)
    # TODO : Not sure what bit value I should choose.
    M = symengine_expr_with_rationals.applyfunc(lambda x: x.n(1000, real=True))

    res = np.empty(sympy_expr.shape[0]*sympy_expr.shape[1], dtype=float)
    for i, val in enumerate(M):
        res[i] = float(val)

    return res.reshape(sympy_expr.shape)


def create_symbol_value_map(constants_name_map, time_varying_name_map):
    """Creates a dictionary mapping the constants and variables to the input
    values specified in BasuMandall2007.

    """

    # These are the Benchmark bicycle [Meijaard, et. al 2007] parameters
    # reexpressed in the Moore 2012 definition.
    bp = bicycle.benchmark_parameters()
    mp = bicycle.benchmark_to_moore(bp)

    # load the input values specified in Table 1 of [BasuMandal2007]_
    basu_input = bicycle.basu_table_one_input()

    # convert the Basu-Mandal values to my coordinates and speeds
    moore_input = bicycle.basu_to_moore_input(basu_input, bp['rR'], bp['lam'])

    # build dictionaries that map the Moore symbolic parameters to the converted
    # Basu-Mandal values
    constant_substitutions = OrderedDict()
    for k, v in mp.items():
        try:
            #exec('constant_substitutions[{}] = v'.format(k))
            #constant_substitutions[sm.Symbol(k)] = v
            constant_substitutions[constants_name_map[k]] = v
        except KeyError:
            print('{} not added to sub dict.'.format(k))
        except NameError:
            print('{} not added to sub dict.'.format(k))

    dynamic_substitutions = {}
    for k, v in moore_input.items():
        try:
            #exec('dynamic_substitutions[{}] = v'.format(k))
            dynamic_substitutions[time_varying_name_map[k]] = v
        except KeyError:
            print('{} not added to sub dict.'.format(k))
        except NameError:
            print('{} not added to sub dict.'.format(k))

    specified_subs = {time_varying_name_map['T4']: 0.0,
                      time_varying_name_map['T6']: 0.0,
                      time_varying_name_map['T7']: 0.0}

    return constant_substitutions, dynamic_substitutions, specified_subs


def compare_cse(expr, args=None):
    """

    Parameters
    ==========
    expr : sympy.Matrix
        I think it has to be a matrix due to some details in the code. Needs to
        be generalized to scalar expression.
    args : List[Symbol,Function(t)]
        If you know all the symbols in the expression already, pass them in.
        Saves some time.

    """

    def find_args(expr):
        args = list(expr.free_symbols)  # get constants and time
        try:  # time not explicit, so remove if there
            args.remove(TIME)
        except ValueError:  # time not present
            pass
        args += me.find_dynamicsymbols(expr)  # get the functions of time
        return list(sm.ordered(args))

    if args is None:
        args = find_args(expr)
    vals = list(np.random.random(len(args)))
    value_dict = {s: v for s, v in zip(args, vals)}

    # TODO : Open an issue in SymPy that tries to lambdify and evalf the
    # bicycle mass matrix, for example.
    # This is where it is getting killed! Can't lambdify M it seems.
    #print('lambdifying the whole expression')
    #eval_expr = sm.lambdify(args, expr)
    #print('Evaluate the whole expression')
    #res = eval_expr(*vals)
    # this also doesn't complete (haven't waited for it to complete)
    #res = sm.matrix2numpy(expr.evalf(subs=value_dict), dtype=float)

    res = sm.matrix2numpy(expr.xreplace(value_dict), dtype=float)

    replacements, reduced_expr = sm.cse(expr)

    for (var, sub_expr) in replacements:
        sub_expr_args = find_args(sub_expr)
        sub_expr_vals = [value_dict[s] for s in sub_expr_args]
        eval_sub_expr = sm.lambdify(sub_expr_args, sub_expr)
        value_dict[var] = eval_sub_expr(*sub_expr_vals)

    red_expr_args = find_args(reduced_expr[0])
    red_expr_vals = [value_dict[s] for s in red_expr_args]
    eval_reduced_expr = sm.lambdify(red_expr_args, reduced_expr[0])
    cse_res = eval_reduced_expr(*red_expr_vals)

    np.testing.assert_allclose(res, cse_res)


def compare_numerically(expr1, expr2, n=10):
    """Compares two SymPy expressions by evaluting with a set of random
    floating point inputs using lambdify."""

    time_varying_symbols1 = me.find_dynamicsymbols(expr1)
    constants1 = expr1.free_symbols
    constants1.remove(TIME)
    args1 = tuple(time_varying_symbols1) + tuple(constants1)

    time_varying_symbols2 = me.find_dynamicsymbols(expr2)
    constants2 = expr2.free_symbols
    constants2.remove(TIME)
    args2 = tuple(time_varying_symbols2) + tuple(constants2)

    # Make the arguments the union of args for both expressions.
    args = tuple(sm.ordered(set(args1 + args2)))

    eval_expr1 = sm.lambdify(args, expr1)
    eval_expr2 = sm.lambdify(args, expr2)

    for i in range(n):
        input_vals = np.random.random(len(args))
        np.testing.assert_allclose(eval_expr1(*input_vals),
                                   eval_expr2(*input_vals))


def decompose_linear_parts(F, *x):
    """Returns the linear coefficient matrices associated with the provided
    vectors and the remainder vector. F must be able to be put into the
    following form:

    F = A1*x1 + A2*x2 + ... + An*xm + B = 0

    - F : n x 1 vector of expressions
    - Ai : n x pi matrix of expressions
    - xi : pi x 1 vector of variables
    - pi : length of vector xi
    - m : number of xi vectors
    - B : n x 1 vector of expressions

    Parameters
    ==========
    F : Matrix, shape(n, 1)
        Column matrix of expressions that linearly depend on entires of
        x1,...,xm.
    x : Sequence[Expr]
        Column matrices representing x1,...,xm.

    Returns
    =======
    Ai, ..., An : Matrix
    B : Matrix, shape(n, 1)

    Notes
    =====
    If xi = xj', then make sure xj'is passed in first to guarantee proper
    replacement.

    """
    F = sm.Matrix(F)
    matrices = []
    for xi in x:
        Ai = F.jacobian(xi)
        matrices.append(Ai)
        repl = {xij: 0 for xij in xi}
        F = F.xreplace(repl)  # remove Ai*xi from F
    matrices.append(F)
    return tuple(matrices)


def formulate_equations_motion(newtonian_frame,
                               bodies,
                               generalized_coordinates, # q
                               generalized_speed_defs,  # u: f(q', q, t)
                               independent_gen_speeds,  # uI
                               # following should be optional
                               dependent_gen_speeds,  # uD
                               motion_constraints,  # G(u, q, t) = 0
                               loads,
                               sub_explicit_gen_dep_speeds=False,
                               nonminimal_form=True):
    """Returns the mass matrix M (coefficients to the generalized
    accelerations) and the forcing vector F (terms that are not functions of
    the generalized accelerations. The result is in one of two forms depending
    if ``nonminimal_form`` is true or not.

    M(q,t)*u'(t) = F(u,q,t)

    of if nonminimal_form

    M(q,t)*uI' = F(u,q,t)

    netownian_frame : ReferenceFrame
        Newtonian reference frame in which the Kane's equations are being
        calculated with respect to.
    bodies : Sequence[Union[Particle, RigidBody]]
        All of the particles and rigid bodies that make up the system S.
    generalized_coordinates : Sequence[Function(t)], len(q)
    generalized_speed_defs : Mapping[Function(t), Expr], len(q)
    independent_gen_speeds : Sequence[Function(t)], len(p)
    dependent_gen_speeds : Sequence[Function(t)], len(m)
    motion_constraints : Sequence[Expr], len(m)
    loads : Mapping[Union[Point, Frame], Vector]
        Mapping of points to their resultant applied force vector and reference
        frames to their resultant applied torque vector.
    sub_explicit_gen_dep_speeds : boolean
    nonminimal_form : boolean

    Returns
    =======
    M : Matrix, shape(n, n) or shape(p, p)
        Matrix containing the coefficients of the generalized accelerations.
    F : Matrix, shape(n, 1) or shape(p, 1)

    Notes
    =====

    Given the p nonohonomic Kane's equations:

    Fr(u,q,t) + Fr*(u',u,q,t) = F(u,q,t) + A_Fs(u,q,t)*u' + B(u, q, t) = 0

    And the m time derivatives of the nonholonomic motion constraints:

    G'(u',u,q,t) = A_G(u,q,t)*u' + B_G(u,q,t) = 0

    the equations can be stacked and written in this form:

    M*u' = F

    where

    M = [A_Fs]
        [A_G]

    F = [-Fr - B]
        [-B_G]

    uI : p independent generalized speeds
    uD : m dependent generalized speeds
    u = [uI]
        [uD]

    """
    t = TIME

    N = newtonian_frame
    q = generalized_coordinates
    uI = independent_gen_speeds
    uD = dependent_gen_speeds
    u = uI + uD

    print('Solving for the time derivatives of the generalized speeds.')
    kin_diff_map = solve_for_qdots(q, generalized_speed_defs)

    print('Solving for the dependent generalized speeds')
    A_GuD, B_G = decompose_linear_parts(motion_constraints, uD)
    uD_of_uI = A_GuD.LUsolve(-B_G)

    # TODO : For some reason if the below lines are used instead of the above
    # two, I get different results in the final ODE evaluation output for the
    # bicycle. But they seem to compare numerically fine.
    #AI, AD, BG = decompose_linear_parts(motion_constraints, uI, uD)
    #uD_of_uI2 = AD.LUsolve(-AI*sm.Matrix(uI) - BG)
    #compare_numerically(uD_of_uI, uD_of_uI2, n=1000)

    if sub_explicit_gen_dep_speeds:
        uD_repl = dict(zip(uD, uD_of_uI))
    else:

        # dependent generalized speeds should be functions of uD(uI, q, t)
        args = tuple(me.find_dynamicsymbols(uD_of_uI))
        print('Partials of the dependent speed expressions with respect to the indepdendent speeds.')
        partials_of_uD = uD_of_uI.jacobian(uI).xreplace(kin_diff_map)
        uD_repl, uD_func_repl, temp_partial_repl, partial_repl = \
            partials_replacements(args, partials_of_uD, uI, uD)

    print('Generating the generalized active forces')
    Fr = generalized_active_forces(N, uI, bodies, loads,
                                   dependent_generalized_speeds=uD_repl)
    print('Generating the generalized inertia forces')
    Frstar = generalized_inertia_forces(N, uI, kin_diff_map, bodies,
                                        dependent_generalized_speeds=uD_repl)

    if not sub_explicit_gen_dep_speeds:
        Fr = Fr.xreplace(temp_partial_repl).xreplace(uD_func_repl)
        Frstar = Frstar.xreplace(temp_partial_repl).xreplace(uD_func_repl)

    print('Formulating the mass matrix form of the equations of motion.')
    if nonminimal_form:
        M, F = formulate_nonmin_mass_matrix_form(Frstar, Fr, u, kin_diff_map,
                                                 motion_constraints)
    else:
        M, F = formulate_mass_matrix_form(Frstar, Fr, uI, uD, kin_diff_map,
                                          motion_constraints)

    if not sub_explicit_gen_dep_speeds:
        M = M.xreplace(partial_repl)
        F = F.xreplace(partial_repl)

    return M, F


def formulate_mass_matrix_form(Frstar, Fr, uI, uD, kin_diff_map,
                               nonholonomic):

    t = TIME

    print('Taking the time derivative of the nonholomic constraints.')
    A_GuI, A_GuD, B_G = decompose_linear_parts(nonholonomic, uI, uD)
    A_GuI_dot = A_GuI.diff(t).xreplace(kin_diff_map)
    A_GuD_dot = A_GuD.diff(t).xreplace(kin_diff_map)
    B_G_dot = B_G.diff(t).xreplace(kin_diff_map)

    print('Decomposing F*')
    A_FstarI, A_FstarD, B_Fstar = decompose_linear_parts(Frstar, uI, uD)

    print('B_F*')
    print(list(sm.ordered(me.find_dynamicsymbols(B_Fstar))))
    print('A_F*I')
    print(list(sm.ordered(me.find_dynamicsymbols(A_FstarI))))
    print('A_F*D')
    print(list(sm.ordered(me.find_dynamicsymbols(A_FstarD))))

    print('Calculate the inverse of A_GuD')
    A_GuD_inv = inv_of_3_by_3(A_GuD)

    uI = sm.Matrix(uI)
    uD = sm.Matrix(uD)

    print("Matrix multiply to get M and F from M*u' = F")
    A = A_FstarD * A_GuD_inv
    M = -A*A_GuI + A_FstarI
    F = A*(A_GuD_dot*uD + A_GuI_dot*uI + B_G_dot) - B_Fstar - Fr

    return M, F


def formulate_nonmin_mass_matrix_form(Frstar, Fr, u, kin_diff_map,
                                      nonholonomic):
    """Returns the coefficient matrix for all generalized accelerations and the
    remaining right hand side of::

       M*u' = F

    where::

       u = [uI]
           [uD]

    Parameters
    ==========
    Frstar : Matrix, shape(p, 1)
    Fr : Matrix, shape(p, 1)
    u : Matrix, shape(n, 1)
    kin_diff_map : Mapping[Derivative(Function(t), t), Expr]
    nonholonomix : Matrix, shape(m, 1)

    Returns
    =======
    M : Matrix, shape(n, n)
    F : Matrix, shape(n, 1)

    """

    u = sm.Matrix(u)
    u_dot = u.diff(TIME)

    G_dot = sm.Matrix(nonholonomic).diff(TIME).xreplace(kin_diff_map)
    A_Gdot, B_Gdot = decompose_linear_parts(G_dot, u_dot)
    A_Fstar, B_Fstar = decompose_linear_parts(Frstar, u_dot)

    MI = A_Fstar  # p x n
    MD = A_Gdot  # m x n
    FI = -B_Fstar - Fr  # p x 1
    FD = -B_Gdot  # m x 1
    M = MI.col_join(MD)
    F = FI.col_join(FD)

    return M, F


def generalized_active_forces(newtonian_frame, generalized_speeds, bodies,
                              loads, dependent_generalized_speeds=None):
    """Returns a column matrix containing p expressions for generalized active
    forces. p is the number of degrees of freedom of the system.

    Parameters
    ==========
    newtonian_frame : ReferenceFrame
        Newtonian reference frame in which the Kane's equations are being
        calculated with respect to.
    generalized_speeds : Sequence[Function], len(p)
        Arbitrary functions of time that represent the independent generalized
        speeds.
    bodies : Iterable[Union[Particle, RigidBody]]
        All of the particles and rigid bodies that make up the system S.
    loads : Mapping[Union[Point, Frame], Vector]
        Mapping of points to their resultant applied force vector and reference
        frames to their resultant applied torque vector.
    dependent_generalized_speeds : Mapping
        Mapping of generalized speed functions of time to either expressions or
        arbitrary functions of the independent speeds.

    Returns
    =======
    Fr : Matrix, shape(p, 1)
        A matrix containing the p expressions for the generalized active
        forces.

    Notes
    =====

    It is expected that the velocities and accelerations for all particles and
    rigid bodies are expressed in terms of the independent generalized speeds.

    """

    Fr = []

    for ur in generalized_speeds:
        rth_generalized_force = sm.S(0)
        for body in bodies:
            contribution = sm.S(0)
            if isinstance(body, me.Particle):
                point = body.point
            if isinstance(body, me.RigidBody):
                point = body.masscenter
                frame = body.frame
                ang_vel = frame.ang_vel_in(newtonian_frame)
                if dependent_generalized_speeds is not None:
                    ang_vel = ang_vel.subs(dependent_generalized_speeds)
                par_ang_vel = ang_vel.diff(ur, newtonian_frame)
                try:
                    torque = loads[frame]
                except KeyError:  # no contribution from this torque
                    pass
                else:
                    contribution += par_ang_vel.dot(torque)
            lin_vel = point.vel(newtonian_frame)
            if dependent_generalized_speeds is not None:
                lin_vel = lin_vel.subs(dependent_generalized_speeds)
            par_lin_vel = lin_vel.diff(ur, newtonian_frame)
            try:
                force = loads[point]
            except KeyError:  # no contribution from this force
                pass
            else:
                contribution += par_lin_vel.dot(force)
            rth_generalized_force += contribution
        Fr.append(rth_generalized_force)

    return sm.Matrix(Fr)


def generalized_inertia_forces(newtonian_frame, generalized_speeds,
                               kinematical_differential_eqs, bodies,
                               dependent_generalized_speeds=None):
    """Returns a column matrix containing p expressions for generalized
    inertial forces. p is the number of degrees of freedom of the system.

    Parameters
    ==========
    newtonian_frame : ReferenceFrame
        Newtonian reference frame in which the Kane's equations are being
        calculated with respect to.
    generalized_speeds : Sequence[Function], len(p)
        Arbitrary functions of time that represent the independent generalized
        speeds.
    kinematical_differential_eqs : Mappping
        A mapping of time derivatives of generalized coordinates to expressions
        that are functions of the generalized coordinates and genearlized
        speeds.
    bodies : Iterable[Union[Particle, RigidBody]]
        All of the particles and rigid bodies that make up the system S.
    dependent_generalized_speeds : Mapping
        Mapping of generalized speed functions of time to either expressions or
        arbitrary functions of the independent speeds.

    Returns
    =======
    Frstar : Matrix, shape(p, 1)
        A matrix containing the p expressions for the generalized inertial
        forces.

    Notes
    =====

    It is expected that the velocities and accelerations for all particles and
    rigid bodies are expressed in terms of the independent generalized speeds.

    """

    Frstar = []

    for ur in generalized_speeds:
        rth_generalized_force = sm.S(0)
        for body in bodies:
            contribution = sm.S(0)
            if isinstance(body, me.Particle):
                point = body.point
            if isinstance(body, me.RigidBody):
                point = body.masscenter
                frame = body.frame
                inertia = body.central_inertia
                ang_vel = frame.ang_vel_in(newtonian_frame)
                if dependent_generalized_speeds is not None:
                    ang_vel = ang_vel.subs(dependent_generalized_speeds)
                ang_acc = frame.ang_acc_in(newtonian_frame)
                # eliminate qdots
                ang_acc = ang_acc.subs(kinematical_differential_eqs)
                par_ang_vel = ang_vel.diff(ur, newtonian_frame)
                inertial_torque = (ang_acc.dot(inertia) +
                                   ang_vel.cross(inertia.dot(ang_vel)))
                contribution += par_ang_vel.dot(-inertial_torque)
            mass = body.mass
            lin_acc = point.acc(newtonian_frame)
            # eliminate qdots
            lin_acc = lin_acc.subs(kinematical_differential_eqs)
            lin_vel = point.vel(newtonian_frame)
            if dependent_generalized_speeds is not None:
                lin_vel = lin_vel.subs(dependent_generalized_speeds)
            par_lin_vel = lin_vel.diff(ur, newtonian_frame)
            contribution += par_lin_vel.dot(-mass*lin_acc)
            rth_generalized_force += contribution
        Frstar.append(rth_generalized_force)

    return sm.Matrix(Frstar)


def inv_of_3_by_3(matrix):
    """Returns the inverse of a 3x3 matrix. The matrix must not be singular.

    Parameters
    ==========
    matrix : Matrix, shape(3, 3)

    Returns
    =======
    invA : Matrix, shape(3, 3)

    """

    a, b, c, d, e, f, g, h, k = matrix

    det = a*(e*k - f*h) - b*(d*k - f*g) + c*(d*h - e*g)

    invA = sm.zeros(3, 3)
    invA[0, 0] = (e*k - f*h) / det
    invA[0, 1] = -(b*k - c*h) / det
    invA[0, 2] = (b*f - c*e) / det
    invA[1, 0] = -(d*k - f*g) / det
    invA[1, 1] = (a*k - c*g) / det
    invA[1, 2] = -(a*f - c*d) / det
    invA[2, 0] = (d*h - e*g) / det
    invA[2, 1] = -(a*h - b*g) / det
    invA[2, 2] = (a*e - b*d) / det

    return invA


def partials_replacements(args, partials_of_uD, uI, uD):
    """Returns several dictionaries useful for managing the dummy functions for
    the partial derivatives of the dependent speeds with respect to the
    independent speeds.

    Parameters
    ==========
    args : Tuple[Union[Function(t), Symbol]]
    partials_of_uD : Matrix, shape(m, p)
    uI : Tuple[Function(t)]
    uD : Tuple[Function(t)]
    """

    uD_funcs = [sm.Function(uDi.name)(*args) for uDi in uD]
    uD_repl = dict(zip(uD, uD_funcs))  # uD(t): uD(uI, q, t)
    uD_func_repl = dict(zip(uD_funcs, uD))  # uD(uI, q, t): uD(t)
    temp_partial_repl = {}
    partial_repl = {}
    for i, uDi in enumerate(uD):
        for j, uIj in enumerate(uI):
            s = me.dynamicsymbols('d{}d{}'.format(uDi.name, uIj.name))
            temp_partial_repl[uD_funcs[i].diff(uIj)] = s
            partial_repl[s] = partials_of_uD[i, j]
    return uD_repl, uD_func_repl, temp_partial_repl, partial_repl


def solve_for_qdots(generalized_coordinates, generalized_speed_definitions):
    """Returns a mapping of generalized coordinate time derivatives to
    expressions of the generalized speeds, coordinates, and time.

    Parameters
    ==========
    generalized_coordinates : Sequence[Function(t)], len(q)
        Arbitrary functions of time representing all of the generalized
        coordinates.
    gen_speed_definitions : Mapping[Function(t), Expr], len(q)
        Mapping of arbitrary functions of time representing the generalized
        speeds to expressions involving the generalized coordinates and their
        time derivatives, i.e. u = f(q', q, t).

    Returns
    =======
    qdots : Mapping[Derivative(Function(t), t), Expr]
        Mapping of the time derivatives of the generalized coordinates to
        expressions involving the generalized speeds, i.e. q' = f(u, q, t).

    """
    t = TIME

    qdot = [qi.diff(t) for qi in generalized_coordinates]

    # the order of the expressions in K does not matter here
    K = sm.Matrix([ui - expr for ui, expr in
                   generalized_speed_definitions.items()])

    repl = {qdoti: 0 for qdoti in qdot}

    # Kinematical differential equations should be linear in the qdot's and the
    # u's. Both of these are true:
    # 1) K = A_Kqd(q,t)*q(t)' + A_Ku(q,t)*u + B_K(q,t) = 0
    # 2) K = A_Kqd(q,t)*q(t)' + B_K(u, q, t) = 0
    # 2) is all that is needed to solve uniquely for the qdot's with this
    # linear system:
    # A_Kqd*q' = -B_K

    A_Kqd = K.jacobian(qdot)
    B_K = K.xreplace(repl)

    qdot_exprs = A_Kqd.LUsolve(-B_K)

    return dict(zip(qdot, qdot_exprs))


def write_matrix_to_file(expr, filename, funcs_of_time=None):

    if funcs_of_time is None:
        funcs_of_time = me.find_dynamicsymbols(expr)
    syms_subs = {f: sm.Symbol(f.name) for f in funcs_of_time}
    expr_with_syms = expr.xreplace(syms_subs)
    expr_srepr = sm.srepr(expr_with_syms)
    with open(filename, 'w') as f:
        f.write(expr_srepr)


def xreplace_and_solve_linear_system(A, b, value_map):
    """Solves Ax=b by substituting values from value_map into A and b and
    solving with NumPy."""

    num_A = sm.matrix2numpy(A.xreplace(value_map), dtype=float)
    num_b = sm.matrix2numpy(b.xreplace(value_map), dtype=float).flatten()

    return np.linalg.solve(num_A, num_b)
