#!/usr/bin/env python

"""This file derives the non-linear equations of motion of the Carvallo-Whipple
bicycle model ([Carvallo1899]_, [Whippl1899]_) following the description and
nomenclature in [Moore2012]_. The resulting equations of motion are compared to
the canonical values presented in [BasuMandal2007]_.

References
==========

.. [Whipple1899] Whipple, Francis J. W. "The Stability of the Motion of a
   Bicycle." Quarterly Journal of Pure and Applied Mathematics 30 (1899):
   312–48.
.. [Carvallo1899] Carvallo, E. Théorie Du Mouvement Du Monocycle et de La
   Bicyclette. Paris, France: Gauthier- Villars, 1899.
.. [Moore2012] Moore, Jason K. "Human Control of a Bicycle." Doctor of
   Philosophy, University of California, 2012.
   http://moorepants.github.io/dissertation.
.. [Meijaard2007] Meijaard, J. P., Jim M. Papadopoulos, Andy Ruina, and A. L.
   Schwab. "Linearized Dynamics Equations for the Balance and Steer of a
   Bicycle: A Benchmark and Review." Proceedings of the Royal Society A:
   Mathematical, Physical and Engineering Sciences 463, no. 2084 (August 8,
   2007): 1955–82.
.. [Basu-Mandal2007] Basu-Mandal, Pradipta, Anindya Chatterjee, and J.M
   Papadopoulos. "Hands-Free Circular Motions of a Benchmark Bicycle."
   Proceedings of the Royal Society A: Mathematical, Physical and Engineering
   Sciences 463, no. 2084 (August 8, 2007): 1983–2003.
   https://doi.org/10.1098/rspa.2007.1849.

"""

import os
import sys

# NOTE : temporary hack so that python pydy/examples/bicycle/whipple.py uses
# the local installed pydy.
sys.path.append(os.path.dirname(__file__))

from pydy.codegen.ode_function_generators import CythonODEFunctionGenerator
import numpy as np
import sympy as sm
import sympy.physics.mechanics as mec

from utils import (
    ReferenceFrame,
    create_basu_output_from_moore_output,
    create_moore_input_from_basu_input,
)

# NOTE : The default cache size is sometimes too low for these large expression
# operations. This potentially helps.
os.environ['SYMPY_CACHE_SIZE'] = '6000'


def setup_symbolics():
    """Returns the fully defined symbolics of the system ready for formulation
    of the equations of motion."""

    ##################
    # Reference Frames
    ##################

    print('Defining reference frames.')

    # Newtonian Frame
    N = ReferenceFrame('N')
    # Yaw Frame
    A = ReferenceFrame('A')
    # Roll Frame
    B = ReferenceFrame('B')
    # Rear Frame
    C = ReferenceFrame('C')
    # Rear Wheel Frame
    D = ReferenceFrame('D')
    # Front Frame
    E = ReferenceFrame('E')
    # Front Wheel Frame
    F = ReferenceFrame('F')

    ####################################
    # Generalized Coordinates and Speeds
    ####################################

    # All the following are a function of time.
    t = mec.dynamicsymbols._t

    # q1: perpendicular distance from the n2> axis to the rear contact
    #     point in the ground plane
    # q2: perpendicular distance from the n1> axis to the rear contact
    #     point in the ground plane
    # q3: frame yaw angle
    # q4: frame roll angle
    # q5: frame pitch angle
    # q6: rear wheel rotation angle
    # q7: steering rotation angle
    # q8: front wheel rotation angle
    # q9: perpendicular distance from the n2> axis to the front contact
    #     point in the ground plane
    # q10: perpendicular distance from the n1> axis to the front contact
    #     point in the ground plane

    print('Defining time varying symbols.')

    q1, q2, q3, q4 = mec.dynamicsymbols('q1 q2 q3 q4')
    q5, q6, q7, q8 = mec.dynamicsymbols('q5 q6 q7 q8')

    u1, u2, u3, u4 = mec.dynamicsymbols('u1 u2 u3 u4')
    u5, u6, u7, u8 = mec.dynamicsymbols('u5 u6 u7 u8')

    # u11: auxiliary speed to determine the rear tire lateral force
    # u12: auxiliary speed to determine the front tire lateral force
    u9, u10, u11, u12 = mec.dynamicsymbols('u9, u10, u11, u12')

    #################################
    # Orientation of Reference Frames
    #################################

    print('Orienting frames.')

    # The following defines a 3-1-2 Tait-Bryan rotation with yaw (q3), roll
    # (q4), pitch (q5) angles to orient the rear frame relative to the ground.
    # The front frame is then rotated through the steer angle (q7) about the
    # rear frame's 3 axis.

    # rear frame yaw
    A.orient(N, 'Axis', (q3, N['3']))
    # rear frame roll
    B.orient(A, 'Axis', (q4, A['1']))
    # rear frame pitch
    C.orient(B, 'Axis', (q5, B['2']))
    # front frame steer
    E.orient(C, 'Axis', (q7, C['3']))

    ###########
    # Constants
    ###########

    print('Defining constants.')

    # geometry
    # rf: radius of front wheel
    # rr: radius of rear wheel
    # d1: the perpendicular distance from the steer axis to the center
    #     of the rear wheel (rear offset)
    # d2: the distance between wheels along the steer axis
    # d3: the perpendicular distance from the steer axis to the center
    #     of the front wheel (fork offset)
    # l1: the distance in the c1> direction from the center of the rear
    #     wheel to the frame center of mass
    # l2: the distance in the c3> direction from the center of the rear
    #     wheel to the frame center of mass
    # l3: the distance in the e1> direction from the front wheel center to
    #     the center of mass of the fork
    # l4: the distance in the e3> direction from the front wheel center to
    #     the center of mass of the fork

    rf, rr = sm.symbols('rf, rr')
    d1, d2, d3 = sm.symbols('d1, d2, d3')
    l1, l2, l3, l4 = sm.symbols('l1, l2, l3, l4')

    # acceleration due to gravity
    g = sm.symbols('g')

    # mass
    mc, md, me, mf = sm.symbols('mc, md, me, mf')

    # inertia components
    ic11, ic22, ic33, ic31 = sm.symbols('ic11, ic22, ic33, ic31')
    id11, id22 = sm.symbols('id11, id22')
    ie11, ie22, ie33, ie31 = sm.symbols('ie11, ie22, ie33, ie31')
    if11, if22 = sm.symbols('if11, if22')

    ###########
    # Specified
    ###########

    # control torques
    # T4 : roll torque
    # T6 : rear wheel torque
    # T7 : steer torque
    # Fr : rear lateral force
    # Ff : front lateral force
    T4, T6, T7, Fr, Ff = mec.dynamicsymbols('T4, T6, T7, Fr, Ff')

    ##################
    # Position Vectors
    ##################

    print('Defining position vectors.')

    # rear wheel contact point
    dn = mec.Point('dn')

    # newtonian origin to rear wheel center
    do = mec.Point('do')
    do.set_pos(dn, -rr*B['3'])

    # rear wheel center to bicycle frame center
    co = mec.Point('co')
    co.set_pos(do, l1*C['1'] + l2*C['3'])

    # rear wheel center to steer axis point
    ce = mec.Point('ce')
    ce.set_pos(do, d1*C['1'])

    # steer axis point to the front wheel center
    fo = mec.Point('fo')
    fo.set_pos(ce, d2*E['3'] + d3*E['1'])

    # front wheel center to front frame center
    eo = mec.Point('eo')
    eo.set_pos(fo, l3*E['1'] + l4*E['3'])

    # front wheel contact point
    fn = mec.Point('fn')
    fn.set_pos(fo, rf*E['2'].cross(A['3']).cross(E['2']).normalize())

    ######################
    # Holonomic Constraint
    ######################

    print('Defining holonomic constraints.')

    # this constraint is enforced so that the front wheel contacts the ground
    holonomic = fn.pos_from(dn).dot(A['3'])

    print('The holonomic constraint is a function of these dynamic variables:')
    print(list(sm.ordered(mec.find_dynamicsymbols(holonomic))))

    ####################################
    # Kinematical Differential Equations
    ####################################

    print('Defining kinematical differential equations.')

    kinematical = [q3.diff(t) - u3,  # yaw
                   q4.diff(t) - u4,  # roll
                   q5.diff(t) - u5,  # pitch
                   q7.diff(t) - u7]  # steer

    ####################
    # Angular Velocities
    ####################

    print('Defining angular velocities.')

    # Note that the wheel angular velocities are defined relative to the frame
    # they are attached to.

    A.set_ang_vel(N, u3*N['3'])  # yaw rate
    B.set_ang_vel(A, u4*A['1'])  # roll rate
    C.set_ang_vel(B, u5*B['2'])  # pitch rate
    D.set_ang_vel(C, u6*C['2'])  # rear wheel rate
    E.set_ang_vel(C, u7*C['3'])  # steer rate
    F.set_ang_vel(E, u8*E['2'])  # front wheel rate

    ###################
    # Linear Velocities
    ###################

    print('Defining linear velocities.')

    # rear wheel contact stays in ground plane and does not slip
    # TODO : Investigate setting to sm.S(0) and 0.
    # include the auxiliary speed corresponding to the later force
    dn.set_vel(N, u11*B['2'])

    # mass centers
    do.v2pt_theory(dn, N, D)
    co.v2pt_theory(do, N, C)
    ce.v2pt_theory(do, N, C)
    fo.v2pt_theory(ce, N, E)
    eo.v2pt_theory(fo, N, E)

    # front wheel contact velocity
    fn.v2pt_theory(fo, N, F)

    # create a front "yaw" frame that is equivalent to the B frame for the rear
    # wheel.
    # G['1'] lies in the ground plane and points in the direction of the wheel
    # contact path E['2'] X N['3'] gives this unit vector.
    # G['2'] lies in the ground plane and points perpendicular to the wheel
    # contact path. N['3'] X G['1'] gives this unit vector.
    g1_hat = E['2'].cross(N['3'])
    g2_hat = N['3'].cross(g1_hat)
    fn.set_vel(N, fn.vel(N) + u12*g2_hat)

    ####################
    # Motion Constraints
    ####################

    print('Defining nonholonomic constraints.')

    nonholonomic = [fn.vel(N).dot(A['1']),
                    fn.vel(N).dot(A['2']),
                    fn.vel(N).dot(A['3'])]

    #########
    # Inertia
    #########

    print('Defining inertia.')

    # NOTE : You cannot define the wheel inertias with respect to their
    # respective frames because the generalized inertia force calcs will fail
    # because there is no direction cosine matrix relating the wheel frames
    # back to the other reference frames so I define them here with respect to
    # the rear and front frames.

    # NOTE : Changing 0.0 to 0 or sm.S(0) changes the floating point errors.

    Ic = mec.inertia(C, ic11, ic22, ic33, 0.0, 0.0, ic31)
    Id = mec.inertia(C, id11, id22, id11, 0.0, 0.0, 0.0)
    Ie = mec.inertia(E, ie11, ie22, ie33, 0.0, 0.0, ie31)
    If = mec.inertia(E, if11, if22, if11, 0.0, 0.0, 0.0)

    ##############
    # Rigid Bodies
    ##############

    print('Defining the rigid bodies.')

    rear_frame = mec.RigidBody('Rear Frame', co, C, mc, (Ic, co))
    rear_wheel = mec.RigidBody('Rear Wheel', do, D, md, (Id, do))
    front_frame = mec.RigidBody('Front Frame', eo, E, me, (Ie, eo))
    front_wheel = mec.RigidBody('Front Wheel', fo, F, mf, (If, fo))

    bodies = [rear_frame, rear_wheel, front_frame, front_wheel]

    ###########################
    # Generalized Active Forces
    ###########################

    print('Defining the forces and torques.')

    # gravity
    Fco = (co, mc*g*A['3'])
    Fdo = (do, md*g*A['3'])
    Feo = (eo, me*g*A['3'])
    Ffo = (fo, mf*g*A['3'])

    # lateral tire forces
    Fdn = (dn, Fr*B['2'])
    Ffn = (fn, Ff*g2_hat)

    # input torques
    Tc = (C, T4*A['1'] - T6*B['2'] - T7*C['3'])
    Td = (D, T6*C['2'])
    Te = (E, T7*C['3'])

    forces = [Fco, Fdo, Feo, Ffo, Fdn, Ffn, Tc, Td, Te]

    # Manually compute the ground contact velocities.
    kindiffdict = sm.solve(kinematical, [q3.diff(t), q4.diff(t), q5.diff(t),
                                         q7.diff(t)], dict=True)[0]
    u1_def = -rr*(u5 + u6)*sm.cos(q3)
    u1p_def = u1_def.diff(t).xreplace(kindiffdict)
    u2_def = -rr*(u5 + u6)*sm.sin(q3)
    u2p_def = u2_def.diff(t).xreplace(kindiffdict)

    ###############################
    # Prep symbolic data for output
    ###############################

    newto = N
    q_ind = (q3, q4, q7)  # yaw, roll, steer
    q_dep = (q5,)  # pitch
    # rear contact 1 dist, rear contact 2 dist, rear wheel angle, front wheel angle
    q_ign = (q1, q2, q6, q8)
    u_ind = (u4, u6, u7)  # roll rate, rear wheel rate, steer rate
    u_dep = (u3, u5, u8)  # yaw rate, pitch rate, front wheel rate
    u_aux = (u11, u12)
    const = (d1, d2, d3, g, ic11, ic22, ic31, ic33, id11, id22, ie11, ie22,
             ie31, ie33, if11, if22, l1, l2, l3, l4, mc, md, me, mf, rf, rr)
    speci = (T4, T6, T7, Fr, Ff)
    holon = [holonomic]
    nonho = tuple(nonholonomic)
    us = tuple(sm.ordered((u1, u2) + u_ind + u_dep + u_aux))
    qs = tuple(sm.ordered(q_ign + q_ind + q_dep))
    # TODO : Reduced these to work with formulate_equations_motion().
    spdef = {ui: qi.diff(t) for ui, qi in zip((u3, u4, u5, u7),
                                              (q3, q4, q5, q7))}
    exdef = {u1: u1_def, u2: u2_def, u1.diff(t): u1p_def, u2.diff(t): u2p_def}

    system_symbolics = {
        'auxiliary speeds': u_aux,
        'bodies': tuple(bodies),
        'constants': const,
        'dependent generalized coordinates': q_dep,
        'dependent generalized speeds': u_dep,
        'extra definitions': exdef,
        'generalized coordinates': qs,
        'generalized speeds': us,
        'holonomic constraints': holon,
        'ignorable coordinates': q_ign,
        'independent generalized coordinates': q_ind,
        'independent generalized speeds': u_ind,
        'kinematical differential equations': kinematical,
        'loads': tuple(forces),
        'newtonian reference frame': newto,
        'nonholonomic constraints': nonho,
        'specified quantities': speci,
        'speed definitions': spdef,
        'time': t,
    }

    return system_symbolics

symbolics = setup_symbolics()

###############
# Kane's Method
###############

print("Generating Kane's equations.")

kane = mec.KanesMethod(
    symbolics['newtonian reference frame'],
    symbolics['independent generalized coordinates'],
    symbolics['independent generalized speeds'],
    kd_eqs=symbolics['kinematical differential equations'],
    q_dependent=symbolics['dependent generalized coordinates'],
    configuration_constraints=symbolics['holonomic constraints'],
    u_dependent=symbolics['dependent generalized speeds'],
    velocity_constraints=symbolics['nonholonomic constraints'],
    u_auxiliary=symbolics['auxiliary speeds'],
)

kane.kanes_equations(symbolics['bodies'], loads=symbolics['loads'])

mass_matrix = kane.mass_matrix
print('The mass matrix is a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(mass_matrix))))

# sub in the kin diffs to eliminate some extraneous derivatives
forcing_vector = kane.forcing.xreplace(kane.kindiffdict())
print('The forcing vector is a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(forcing_vector))))

print('The auxiliary equations are a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(kane.auxiliary_eqs))))

u_aux = sm.Matrix(symbolics['auxiliary speeds'])
lat_forces = sm.Matrix(symbolics['specified quantities'][-2:])
zero_lat_forces = {f: 0 for f in lat_forces}

"""
Should be linear in the forces?

A(q, t)*[Fr] - b(u', u, q, t) = 0
        [Ff]

"""

A = kane.auxiliary_eqs.jacobian(lat_forces)
b = -kane.auxiliary_eqs.xreplace(zero_lat_forces)

print(list(sm.ordered(mec.find_dynamicsymbols(A))))
print(list(sm.ordered(mec.find_dynamicsymbols(b))))

from pydy.codegen.octave_code import OctaveMatrixGenerator

us = symbolics['generalized speeds']
u_dots = [mec.dynamicsymbols(ui.name + 'p') for ui in us]
u_dot_subs = {ui.diff(): upi for ui, upi in zip(us, u_dots)}

b = b.xreplace(u_dot_subs)

xs = list(sm.ordered(mec.find_dynamicsymbols(A).union(mec.find_dynamicsymbols(b))))
ps = list(sm.ordered(A.free_symbols.union(b.free_symbols)))
ps.remove(symbolics['time'])

gen = OctaveMatrixGenerator([xs, ps], [A, b])
gen.write('eval_lat_forces')
