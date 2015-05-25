#!usr/bin/env python

"""This file derives the non-linear equations of motion of the Whipple
bicycle model following the description and nomenclature in [Moore2012]_.

It option depends on DynamicistToolKit to compare to the canonical values of
this problem.

.. [Moore2012] Moore, Jason K. "Human Control of a Bicycle." Doctor of
   Philosophy, University of California, Davis, 2012.
   http://moorepants.github.com/dissertation

"""

# TODO : Make dtk optional.

from collections import OrderedDict

import numpy as np
import sympy as sm
import sympy.physics.mechanics as mec
from pydy.codegen.ode_function_generators import CythonODEFunctionGenerator
from dtk import bicycle

mec.Vector.simp = False

##################
# Reference Frames
##################

print('Defining reference frames.')


class ReferenceFrame(mec.ReferenceFrame):
    """Subclass that enforces the desired unit vector indice style."""

    def __init__(self, *args, **kwargs):

        kwargs.pop('indices', None)
        kwargs.pop('latexs', None)

        lab = args[0].lower()
        tex = '\hat{{{}}}_{}'

        super(ReferenceFrame, self).__init__(*args,
                                             indices=('1', '2', '3'),
                                             latexs=(tex.format(lab, '1'),
                                                     tex.format(lab, '2'),
                                                     tex.format(lab, '3')),
                                             **kwargs)

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
# q6: front wheel rotation angle
# q7: steering rotation angle
# q8: rear wheel rotation angle
# q9: perpendicular distance from the n2> axis to the front contact
#     point in the ground plane
# q10: perpendicular distance from the n1> axis to the front contact
#     point in the ground plane

print('Defining time varying symbols.')

q1, q2, q3, q4 = mec.dynamicsymbols('q1 q2 q3 q4')
q5, q6, q7, q8 = mec.dynamicsymbols('q5 q6 q7 q8')

u1, u2, u3, u4 = mec.dynamicsymbols('u1 u2 u3 u4')
u5, u6, u7, u8 = mec.dynamicsymbols('u5 u6 u7 u8')

#################################
# Orientation of Reference Frames
#################################

print('Orienting frames.')

# TODO : Report this as a SymPy issue. The following fails:
#A.orient(N, 'Axis', (q3, A['3']))

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

rf, rr = sm.symbols('rf rr', real=True, positive=True)
d1, d2, d3 = sm.symbols('d1 d2 d3', real=True)
l1, l2, l3, l4 = sm.symbols('l1 l2 l3 l4', real=True)

# acceleration due to gravity
g = sm.symbols('g', real=True)

# mass
mc, md, me, mf = sm.symbols('mc md me mf', real=True, positive=True)

# inertia
ic11, ic22, ic33, ic31 = sm.symbols('ic11 ic22 ic33 ic31', real=True)
id11, id22 = sm.symbols('id11 id22', real=True)
ie11, ie22, ie33, ie31 = sm.symbols('ie11 ie22 ie33 ie31', real=True)
if11, if22 = sm.symbols('if11 if22', real=True)

###########
# Specified
###########

# control torques
# T4 : roll torque
# T6 : rear wheel torque
# T7 : steer torque
T4, T6, T7 = mec.dynamicsymbols('T4 T6 T7')

##################
# Position Vectors
##################

print('Defining position vectors.')

# newtonian origin
no = mec.Point('no')

# newtonian origin to rear wheel center
do = mec.Point('do')
do.set_pos(no, -rr * B['3'])

# rear wheel center to bicycle frame center
co = mec.Point('co')
co.set_pos(do, l1 * C['1'] + l2 * C['3'])

# rear wheel center to steer axis point
ce = mec.Point('ce')
ce.set_pos(do, d1 * C['1'])

# steer axis point to the front wheel center
fo = mec.Point('fo')
fo.set_pos(ce, d2 * E['3'] + d3 * E['1'])

# front wheel center to front frame center
eo = mec.Point('eo')
eo.set_pos(fo, l3 * E['1'] + l4 * E['3'])

# locate the points fixed on the wheel which instaneously touch the ground
# rear
dn = mec.Point('dn')
dn.set_pos(do, rr * B['3'])
# front
fn = mec.Point('fn')
fn.set_pos(fo, rf * E['2'].cross(A['3']).cross(E['2']).normalize())

######################
# Holonomic Constraint
######################

print('Defining holonomic constraints.')

# this constraint is enforced so that the front wheel contacts the ground
holonomic = fn.pos_from(dn).dot(A['3'])
#holonomic = sm.trigsimp(holonomic)

print('The holonomic constraint is a function of these dynamice variables:')
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

A.set_ang_vel(N, u3 * N['3'])  # yaw rate
B.set_ang_vel(A, u4 * A['1'])  # roll rate
C.set_ang_vel(B, u5 * B['2'])  # pitch rate
D.set_ang_vel(C, u6 * C['2'])  # rear wheel rate
E.set_ang_vel(C, u7 * C['3'])  # steer rate
F.set_ang_vel(E, u8 * E['2'])  # front wheel rate

###################
# Linear Velocities
###################

print('Defining linear velocities.')

# origin is fixed
no.set_vel(N, 0.0 * N['1'])

# mass centers
# THIS CHANGE WAS THE FUCKING MAIN ERROR!
#do.set_vel(N, do.pos_from(no).dt(N))
do.v2pt_theory(no, N, D)
co.v2pt_theory(do, N, C)
ce.v2pt_theory(do, N, C)
fo.v2pt_theory(ce, N, E)
eo.v2pt_theory(fo, N, E)

# wheel contact velocities
dn.set_vel(N, 0.0 * N['1'])
fn.v2pt_theory(fo, N, F)

####################
# Motion Constraints
####################

print('Defining nonholonomic constraints.')

nonholonomic = [fn.vel(N).dot(A['1']),
                fn.vel(N).dot(A['2']),
                fn.vel(N).dot(A['3'])]
# The following is pretty slow.
#nonholonomic = [sm.trigsimp(expr) for expr in nonholonomic]

# TODO : simplify(nh1-nh2) doesn't give me zero
#nh1 = fn.vel(N).dot(A['3'])
#nh2 = holonomic.diff(t).subs(sm.solve(kinematical, [q3.diff(t), q4.diff(t),
                                                    #q5.diff(t), q7.diff(t)],
                                      #dict=True)[0])

print('The nonholonomic constraints are a function of these dynamice variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(sm.Matrix(nonholonomic)))))

#########
# Inertia
#########

print('Defining inertia.')

# NOTE : You cannot define the wheel inertias with respect to their
# respective frames because the generalized inertia force calcs will fail
# because there is no direction cosine matrix relating the wheel frames back
# to the other reference frames so I define them here with respect to the
# rear and front frames.

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
Fco = (co, mc * g * A['3'])
Fdo = (do, md * g * A['3'])
Feo = (eo, me * g * A['3'])
Ffo = (fo, mf * g * A['3'])

# input torques
Tc = (C, T4 * A['1'] - T6 * B['2'] - T7 * C['3'])
Td = (D, T6 * C['2'])
Te = (E, T7 * C['3'])

forces = [Fco, Fdo, Feo, Ffo, Tc, Td, Te]

###############
# Kane's Method
###############

print("Generating Kane's equations.")

# var   : J , G
# yaw   : q3, q1
# roll  : q4, q2
# pitch : q5, q4
# steer : q7, q5
# rwh   : u6, u3
# fwh   : u8, u6

kane = mec.KanesMethod(N,
                       [q3, q4, q7],  # yaw, roll, steer
                       [u4, u6, u7],  # roll rate, rear wheel rate, steer rate
                       kd_eqs=kinematical,
                       q_dependent=[q5],  # pitch angle
                       configuration_constraints=[holonomic],
                       u_dependent=[u3, u5, u8],  # yaw rate, pitch rate, front wheel rate
                       velocity_constraints=nonholonomic)

fr, frstar = kane.kanes_equations(forces, bodies)

mass_matrix = kane.mass_matrix
print('The mass matrix is a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(mass_matrix))))

# sub in the kin diffs to eliminate some extraneous derivatives
forcing_vector = mec.msubs(kane.forcing, kane.kindiffdict())
print('The forcing vector is a function of these dynamic variables:')
print(list(sm.ordered(mec.find_dynamicsymbols(forcing_vector))))

# NOTE : doing this substitution must make the expressions super long
# because it causes all kinds of slowdowns, for example cse is slow.
# u3, u5, u8 are all in forcing_vector and need to be replaced by functions
# of u4, u6, u7
# These are the equations for u3, u5, u8 in terms of the independent speeds:
#u_dep = kane._Ars * kane.u[:3, :]
#print('The dependent speeds are a function of these dynamic variables:')
#print(list(sm.ordered(mec.find_dynamicsymbols(u_dep))))
#
#forcing_vector = mec.msubs(forcing_vector, dict(zip([u3, u5, u8], u_dep)))
# This takes forever.
#print(list(sm.ordered(mec.find_dynamicsymbols(forcing_vector))))

# Validation of non-linear equations

print('Loading numerical input parameters.')

bp = bicycle.benchmark_parameters()
mp = bicycle.benchmark_to_moore(bp)

# load the input values specified in Table 1 of Basu-Mandal2007
basu_input = bicycle.basu_table_one_input()

# convert the Basu-Mandall values to my coordinates and speeds
moore_input = bicycle.basu_to_moore_input(basu_input, bp['rR'], bp['lam'])

constant_substitutions = OrderedDict()
for k, v in mp.items():
    try:
        exec('constant_substitutions[{}] = v'.format(k))
        #constant_substitutions[sm.Symbol(k)] = v
    except NameError:
        print('{} not added to sub dict.'.format(k))

dynamic_substitutions = {}
for k, v in moore_input.items():
    try:
        exec('dynamic_substitutions[{}] = v'.format(k))
    except NameError:
        print('{} not added to sub dict.'.format(k))

specified_subs = {T4: 0.0, T6: 0.0, T7: 0.0}
substitutions = specified_subs.copy()
substitutions.update(constant_substitutions)
substitutions.update(dynamic_substitutions)

# Try substituting values in through SymPy
print('Substituting numerical parameters.')
num_mass_matrix = mec.msubs(mass_matrix, substitutions)
num_forcing_vector = mec.msubs(forcing_vector, substitutions)

print("Solving for x'.")
xd = num_mass_matrix.LUsolve(num_forcing_vector)

print(xd)

# BUGS to report:
# 1. find_dynamicsymbols should deal with Vectors and lists of
# exprs/vectors/etc.

# TODO : cse takes forever with this problem if the udeps are substituted
# into forcing_vector, but may be crucial for speed.
# TODO : The above xd computes but I get singular matrix error for this.
rhs_of_kin_diffs = sm.Matrix([kane.kindiffdict()[k] for k in kane.q.diff(t)])
g = CythonODEFunctionGenerator(forcing_vector,
                               kane.q[:], # q3, q4, q7
                               kane.u[:], # u4, u6, u7
                               constant_substitutions.keys(),
                               mass_matrix=mass_matrix,
                               coordinate_derivatives=rhs_of_kin_diffs,
                               specifieds=[T4, T6, T7],
                               constants_arg_type='array',
                               specifieds_arg_type='array')
print('Generating rhs')
#rhs = g.generate()

state_vals = []
for d in kane.q[:] + kane.u[:]:
    sym_str = str(d)[:-3]
    state_vals.append(moore_input[sym_str])
state_vals = np.array(state_vals)

#xd = rhs(state_vals, 0.0, np.array(constant_substitutions.values()),
         #np.array([0.0, 0.0, 0.0]))
#print(xd)

print("Generating output dictionary.")
# convert the outputs from my model to the Basu-Mandal coordinates
# TODO : raise an issue about not knowing which order the x vector is in with
# reference to M * x' = F
speeds = kane.u[:]
speed_deriv_names = [str(speed)[:-3] + 'p' for speed in speeds]
moore_output = {k: v for k, v in zip(speed_deriv_names, list(xd))}
u1 = -rr * u6 * sm.cos(q3)
u1p = u1.diff(t)
u2 = -rr * u6 * sm.sin(q3)
u2p = u2.diff(t)
moore_output['u1p'] = u1p.subs({u6.diff(t): moore_output['u6p']}).subs(kane.kindiffdict()).subs(substitutions)
moore_output['u2p'] = u2p.subs({u6.diff(t): moore_output['u6p']}).subs(kane.kindiffdict()).subs(substitutions)
moore_output.update(moore_input)

moore_output_basu = bicycle.moore_to_basu(moore_output, bp['rR'], bp['lam'])
basu_output = bicycle.basu_table_one_output()

print("Assertions.")
from numpy import testing

for k, bv in basu_output.items():
    mv = float(moore_output_basu[k])
    try:
        testing.assert_allclose(bv, mv)
    except AssertionError:
        print('Failed: {} is supposed to be {:1.16f} but is {:1.16f}'.format(k, bv, mv))
