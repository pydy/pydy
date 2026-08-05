#!/usr/bin/env python

"""This file derives the non-linear equations of motion of the Carvallo-Whipple
bicycle model ([Carvallo1899]_, [Whippl1899]_) following the description and
nomenclature in [Moore2012]_ and produces Octave functions that calculate the
lateral wheel-ground constraint force for each wheel given the essential
kinematics of the vehicle.

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

"""

import os

import sympy as sm
import sympy.physics.mechanics as mec
from pydy.codegen.octave_code import OctaveMatrixGenerator

# NOTE : The default cache size is sometimes too low for these large expression
# operations. This potentially helps.
os.environ['SYMPY_CACHE_SIZE'] = '6000'


class ReferenceFrame(mec.ReferenceFrame):
    """Subclass that enforces the desired unit vector index style."""

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

print('Defining time varying symbols.')

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
q1, q2, q3, q4 = mec.dynamicsymbols('q1, q2, q3, q4')
q5, q6, q7, q8 = mec.dynamicsymbols('q5, q6, q7, q8')
q10, q11 = mec.dynamicsymbols('q10, q11')

# u1: speed of the rear wheel contact point in the n1> direction
# u2: speed of the rear wheel contact point in the n2> direction
# u3: frame yaw angular rate
# u4: frame roll angular rate
# u5: frame pitch angular rate
# u6: rear wheel rotation angular rate
# u7: steering rotation angular rate
# u8: front wheel rotation angular rate
u1, u2, u3, u4 = mec.dynamicsymbols('u1, u2, u3, u4')
u5, u6, u7, u8 = mec.dynamicsymbols('u5, u6, u7, u8')

# u9: speed of the front wheel contact point in the n1> direction
# u10: speed of the front wheel contact point in the n2> direction
# u11: auxiliary speed to determine the rear tire lateral force
# u12: auxiliary speed to determine the front tire lateral force
u9, u10, u11, u12 = mec.dynamicsymbols('u9, u10, u11, u12')

#################################
# Orientation of Reference Frames
#################################

print('Orienting frames.')

# The following defines a 3-1-2 Tait-Bryan rotation with yaw (q3), roll
# (q4), pitch (q5) angles to orient the rear frame relative to the ground
# (Newtonian frame). The front frame is then rotated through the steer
# angle (q7) about the rear frame's 3 axis. The wheels are not oriented, as
# q6 and q8 end up being ignorable coordinates.

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

# mass for each rigid body: C, D, E, F
mc, md, me, mf = sm.symbols('mc, md, me, mf')

# inertia components for each rigid body: C, D, E, F
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
# Fr : rear wheel-ground contact lateral force
# Ff : front wheel-ground contact lateral force
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

# rear wheel contact stays in ground plane and does not slip but the auxiliary
# speed, u11, is included which corresponds to the later force
dn.set_vel(N, 0)
dn_ = mec.Point('dn')
dn_.set_pos(dn, 0)
dn_.set_vel(N, dn.vel(N) + u11*A['2'])

# mass centers
do.v2pt_theory(dn_, N, D)  # ensures u11 in present in velocities
co.v2pt_theory(do, N, C)
ce.v2pt_theory(do, N, C)
fo.v2pt_theory(ce, N, E)
eo.v2pt_theory(fo, N, E)

# front wheel contact velocity
fn.v2pt_theory(fo, N, F)

# create a front "yaw" frame that is equivalent to the A frame for the rear
# wheel.
# G['1'] lies in the ground plane and points in the direction of the wheel
# contact path E['2'] X A['3'] gives this unit vector.
# G['2'] lies in the ground plane and points perpendicular to the wheel
# contact path. A['3'] X G['1'] gives this unit vector.
g1_hat = E['2'].cross(A['3'])
g2_hat = A['3'].cross(g1_hat)
fn_ = mec.Point('fn')
fn_.set_pos(fn, 0)
fn_.set_vel(N, fn.vel(N) + u12*g2_hat)  # includes u11 and u12

####################
# Motion Constraints
####################

print('Defining nonholonomic constraints.')

nonholonomic = [fn_.vel(N).dot(A['1']),
                fn_.vel(N).dot(A['2']),
                fn_.vel(N).dot(A['3'])]

#########
# Inertia
#########

print('Defining inertia.')

Ic = mec.inertia(C, ic11, ic22, ic33, 0, 0, ic31)
Id = mec.inertia(C, id11, id22, id11, 0, 0, 0)
Ie = mec.inertia(E, ie11, ie22, ie33, 0, 0, ie31)
If = mec.inertia(E, if11, if22, if11, 0, 0, 0)

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

# lateral tire forces, need equal and opposite forces
Fdn = (dn, Fr*A['2'])
Fdn_ = (dn_, -Fr*A['2'])
Ffn = (fn, Ff*g2_hat)
Ffn_ = (fn_, -Ff*g2_hat)

# input torques
Tc = (C, T4*A['1'] - T6*B['2'] - T7*C['3'])
Td = (D, T6*C['2'])
Te = (E, T7*C['3'])

forces = [Fco, Fdo, Feo, Ffo, Fdn, Ffn, Fdn_, Ffn_, Tc, Td, Te]

# Manually compute the ground contact velocities.
kindiffdict = sm.solve(kinematical, [q3.diff(t), q4.diff(t), q5.diff(t),
                                     q7.diff(t)], dict=True)[0]
u1_def = -rr*(u5 + u6)*sm.cos(q3)
u1p_def = u1_def.diff(t).xreplace(kindiffdict)
u2_def = -rr*(u5 + u6)*sm.sin(q3)
u2p_def = u2_def.diff(t).xreplace(kindiffdict)

####################
# Prep symbolic data
####################

newto = N
q_ind = (q3, q4, q7)  # yaw, roll, steer
q_dep = (q5,)  # pitch
# rear contact 1 dist, rear contact 2 dist, rear wheel angle, front wheel angle
q_ign = (q1, q2, q6, q8)
qs = tuple(sm.ordered(q_ign + q_ind + q_dep))

u_ind = (u4, u6, u7)  # roll rate, rear wheel rate, steer rate
u_dep = (u3, u5, u8)  # yaw rate, pitch rate, front wheel rate
u_aux = (u11, u12)
us = tuple(sm.ordered((u1, u2) + u_ind + u_dep + u_aux))

const = (d1, d2, d3, g, ic11, ic22, ic31, ic33, id11, id22, ie11, ie22, ie31,
         ie33, if11, if22, l1, l2, l3, l4, mc, md, me, mf, rf, rr)
speci = (T4, T6, T7, Fr, Ff)
holon = [holonomic]
nonho = tuple(nonholonomic)

###############
# Kane's Method
###############

print("Generating Kane's equations.")

kane = mec.KanesMethod(
    newto,
    q_ind,
    u_ind,
    kd_eqs=kinematical,
    q_dependent=q_dep,
    configuration_constraints=holon,
    u_dependent=u_dep,
    velocity_constraints=nonho,
    u_auxiliary=u_aux,
)

kane.kanes_equations(bodies, loads=forces)

###############
# IMU Equations
###############

u4p, u6p, u7p = mec.dynamicsymbols('u4p, u6p, u7p')
u3p, u5p, u8p = mec.dynamicsymbols('u3p, u5p, u8p')
u_dots = [mec.dynamicsymbols(ui.name + 'p') for ui in us]
u_dot_subs = {ui.diff(): upi for ui, upi in zip(us, u_dots)}

diff_subs = {u3.diff(): u3p,
             u4.diff(): u4p,
             u5.diff(): u5p,
             u6.diff(): u6p,
             u7.diff(): u7p,
             u11.diff(): 0,  # fictitious
             u12.diff(): 0}  # fictitious
non_diff_subs = {u11: 0,  # fictitious
                 u12: 0}  # fictitious

bx, by, bz, ex, ey, ez = sm.symbols('bx, by, bz, ex, ey, ez')

lam = sm.symbols('lambda')
C_lam = C.orientnew('C_lam', 'Axis', (lam, -B['2']))

# point in the rear frame measured from rear wheel center
P = do.locatenew('P', bx*C_lam['1'] + by*C_lam['2'] + bz*C_lam['3'])
P.v2pt_theory(do, N, C_lam)
P.a2pt_theory(do, N, C_lam)

# point in the steered frame measured from front wheel center
Q = fo.locatenew('Q', ex*E['1'] + ey*E['2'] + ez*E['3'])
Q.v2pt_theory(fo, N, E)
Q.a2pt_theory(fo, N, E)

# NOTE : This assumes that the smartphones XYZ axes are aligned with the B and
# E axes when in the no roll, pitch, steer configuration.

eqs = sm.Matrix([
    # angular velocity components of the rear frame expressed in the roll frame
    # because IMU is aligned with the roll frame coordinates
    C.ang_vel_in(N).dot(C_lam['1']),
    C.ang_vel_in(N).dot(C_lam['2']),
    C.ang_vel_in(N).dot(C_lam['3']),
    # body fixed linear acceleration of point P (including gravity) expressed
    # in the roll frame because IMU is aligned with the roll frame coordinates
    (P.acc(N) - g*A['3']).dot(C_lam['1']),
    (P.acc(N) - g*A['3']).dot(C_lam['2']),
    (P.acc(N) - g*A['3']).dot(C_lam['3']),
    # body fixed angular velocity components of the steer frame
    E.ang_vel_in(N).dot(E['1']),
    E.ang_vel_in(N).dot(E['2']),
    E.ang_vel_in(N).dot(E['3']),
    # body fixed linear acceleration components of a point on the steer frame
    (Q.acc(N) - g*A['3']).dot(E['1']),
    (Q.acc(N) - g*A['3']).dot(E['2']),
    (Q.acc(N) - g*A['3']).dot(E['3']),
])

eqs = eqs.xreplace(kindiffdict).xreplace(diff_subs).xreplace(non_diff_subs)

gen = OctaveMatrixGenerator([[q4, q5, q7],
                             [u3, u4, u5, u6, u7],
                             [u3p, u4p, u5p, u6p, u7p],
                             [bx, by, bz, d1, d2, d3, ex, ey, ez, g, lam, rr]],
                            [eqs])
gen.write('eval_imu', path=os.path.dirname(__file__))


###########################
# Generate Octave Functions
###########################

gen = OctaveMatrixGenerator([[q4, q5, q7],
                             [d1, d2, d3, rf, rr]],
                            [sm.Matrix([holonomic])])
gen.write('eval_holonomic', path=os.path.dirname(__file__))

# Create matrices for solving for the dependent speeds.
nonholonomic = sm.Matrix(nonholonomic).xreplace({u11: 0, u12: 0})

print('The nonholonomic constraints a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(nonholonomic))))

A_nh, B_nh = decompose_linear_parts(nonholonomic, [u3, u5, u8])
gen = OctaveMatrixGenerator([[q4, q5, q7],
                             [u4, u6, u7],
                             [d1, d2, d3, rf, rr]],
                            [A_nh, -B_nh])
gen.write('eval_dep_speeds', path=os.path.dirname(__file__))

# Create function for solving for the derivatives of the dependent speeds.
nonholonomic_dot = sm.Matrix(nonholonomic).diff(t).xreplace(kane.kindiffdict())

nonholonomic_dot = nonholonomic_dot.xreplace(u_dot_subs)

print('The derivative of the nonholonomic constraints a function of these '
      'dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(nonholonomic_dot))))

A_pnh, B_pnh = decompose_linear_parts(nonholonomic_dot,
                                      [u3p, u5p, u8p])
gen = OctaveMatrixGenerator([[q4, q5, q7],
                             [u3, u4, u5, u6, u7, u8],
                             [u4p, u6p, u7p],
                             [d1, d2, d3, rf, rr]],
                            [A_pnh, -B_pnh])
gen.write('eval_dep_speeds_derivs', path=os.path.dirname(__file__))

# Create function for solving for the lateral forces.
"""
Should be linear in the forces? Or even always F1 + F2 + ... = 0, i.e.
coefficient is 1?

A(q, t)*[Ff] - b(u', u, q, t) = 0
        [Fr]

"""

aux_eqs = kane.auxiliary_eqs.xreplace({u11: 0, u12: 0}).xreplace(
    u_dot_subs).xreplace(kane.kindiffdict())
print('The auxiliary equations are a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(aux_eqs))))

# TODO: Ff is only in aux_eq[1], so decompose fails when trying to take the
# Jacobian wrt Fr in decompose_lienar_parts. Oddly it doesn't just return 0 for
# that component.
a11 = aux_eqs[0].diff(Ff)
a12 = aux_eqs[0].diff(Fr)
a21 = aux_eqs[1].diff(Ff)
a22 = aux_eqs[1].diff(Fr)
A = sm.Matrix([[a11, a12], [a21, a22]])
b = -aux_eqs.xreplace({Ff: 0, Fr: 0})

print('A is a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(A))))
print('b is a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(b))))

gen = OctaveMatrixGenerator([[q4, q5, q7],
                             [u3, u4, u5, u6, u7, u8],
                             [u3p, u4p, u5p, u6p, u7p, u8p],
                             list(const)],
                            [A, b])
gen.write('eval_lat_forces', path=os.path.dirname(__file__))
