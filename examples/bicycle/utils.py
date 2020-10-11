#!/usr/bin/env python

from collections import OrderedDict

import numpy as np
import sympy as sm
import sympy.physics.mechanics as me
from dtk import bicycle

TIME = me.dynamicsymbols._t


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


def compare_numerically(expr1, expr2, n=10):
    """Compares two SymPy expressions by evaluting with a set of random
    floating point inputs."""

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
    """Returns the mass matrix M (coefficients to the generalized accelerations)
    and the forcing vector F (terms that are not functions of the generalized
    accelerations. The result is in one of two forms depending if
    ``nonminimal_form`` is true or not.

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


def solve_linear_system_with_sympy_subs(A, b, value_map):
    """Solves Ax=b by substituting values from value_map into A and b and
    solving with NumPy."""

    num_A = sm.matrix2numpy(A.xreplace(value_map), dtype=float)
    num_b = sm.matrix2numpy(b.xreplace(value_map), dtype=float).flatten()

    return np.linalg.solve(num_A, num_b)
