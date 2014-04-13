#!usr/bin/env python

# This file derives the non-linear equations of motion of the Whipple bicycle
# model.

import sympy as sym
import sympy.physics.mechanics as mec

# debugging
try:
    from IPython.core.debugger import Tracer
except ImportError:
    pass
else:
    set_trace = Tracer()

mec.Vector.simp = False

##################
# Reference Frames
##################

print('Defining reference frames.')

# Newtonian Frame
N = mec.ReferenceFrame('N', indices=('1', '2', '3'),
        latexs=('\hat{n}_1', '\hat{n}_2', '\hat{n}_3'))
# Yaw Frame
A = mec.ReferenceFrame('A', indices=('1', '2', '3'),
        latexs=('\hat{a}_1', '\hat{a}_2', '\hat{a}_3'))
# Roll Frame
B = mec.ReferenceFrame('B', indices=('1', '2', '3'),
        latexs=('\hat{b}_1', '\hat{b}_2', '\hat{b}_3'))
# Rear Frame
C = mec.ReferenceFrame('C', indices=('1', '2', '3'),
        latexs=('\hat{c}_1', '\hat{c}_2', '\hat{c}_3'))
# Rear Wheel Frame
D = mec.ReferenceFrame('D', indices=('1', '2', '3'),
        latexs=('\hat{d}_1', '\hat{d}_2', '\hat{d}_3'))
# Front Frame
E = mec.ReferenceFrame('E', indices=('1', '2', '3'),
        latexs=('\hat{e}_1', '\hat{e}_2', '\hat{e}_3'))
# Front Wheel Frame
F = mec.ReferenceFrame('F', indices=('1', '2', '3'),
        latexs=('\hat{f}_1', '\hat{f}_2', '\hat{f}_3'))

####################################
# Generalized Coordinates and Speeds
####################################

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
q1d, q2d, q3d, q4d = mec.dynamicsymbols('q1 q2 q3 q4', 1)
q5d, q6d, q7d, q8d = mec.dynamicsymbols('q5 q6 q7 q8', 1)

u1, u2, u3, u4 = mec.dynamicsymbols('u1 u2 u3 u4')
u5, u6, u7, u8 = mec.dynamicsymbols('u5 u6 u7 u8')
u1d, u2d, u3d, u4d = mec.dynamicsymbols('u1 u2 u3 u4', 1)
u5d, u6d, u7d, u8d = mec.dynamicsymbols('u5 u6 u7 u8', 1)

ua1, ua2, ua3, ua4, ua5, ua6 = mec.dynamicsymbols('ua1  ua2  ua3  ua4  ua5 ua6')

# TODO : switch to some more compact code
#q = mec.dynamicsymbols('q:8')
#qd = mec.dynamicsymbols('q:8', 1)
#u = mec.dynamicsymbols('u:8')
#ud = mec.dynamicsymbols('u:8', 1)
#ua = mec.dynamicsymbols('ua:6')

#################################
# Orientation of Reference Frames
#################################

print('Orienting frames.')

# TODO : report this issue, or fix it. The following fails.
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

rf, rr = sym.symbols('rf rr')
d1, d2, d3 = sym.symbols('d1 d2 d3')
l1, l2, l3, l4 = sym.symbols('l1 l2 l3 l4')

# acceleration due to gravity
g = sym.symbols('g')

# mass
mc, md, me, mf = sym.symbols('mc md me mf')

# inertia
ic11, ic22, ic33, ic31 = sym.symbols('ic11 ic22 ic33 ic31')
id11, id22 = sym.symbols('id11 id22')
ie11, ie22, ie33, ie31 = sym.symbols('ie11 ie22 ie33 ie31')
if11, if22 = sym.symbols('if11 if22')

# time
t = sym.symbols('t')

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

####################################
# Kinematical Differential Equations
####################################

print('Defining kinematical differential equations.')

kinematical = [q3d - u3, # yaw
               q4d - u4, # roll
               q5d - u5, # pitch
               q7d - u7] # steer

####################
# Angular Velocities
####################

print('Defining angular velocities.')

A.set_ang_vel(N, u3 * N['3']) # yaw rate
B.set_ang_vel(A, u4 * A['1']) # roll rate
C.set_ang_vel(B, u5 * B['2']) # pitch rate
D.set_ang_vel(C, u6 * C['2']) # rear wheel rate
E.set_ang_vel(C, u7 * C['3']) # steer rate
F.set_ang_vel(E, u8 * E['2']) # front wheel rate

###################
# Linear Velocities
###################

print('Defining linear velocities.')

# origin is fixed
no.set_vel(N, 0.0 * N['1'])

# mass centers
do.set_vel(N, do.pos_from(no).dt(N))
co.v2pt_theory(do, N, C)
ce.v2pt_theory(do, N, C)
fo.v2pt_theory(ce, N, E)
eo.v2pt_theory(fo, N, E)

# wheel contact velocities
dn.set_vel(N, 0)
fn.v2pt_theory(fo, N, F)

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

# note that you cannot define the wheel inertia's with respect to their
# respective frames because the generalized inertia force calcs will fail
# because there is no direction cosine matrix relating the wheel frames back to
# the other reference frames so I define them here with respect to the rear and
# front frames

print('Defining inertia.')

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

kane = mec.KanesMethod(N,
                       [q3, q4, q7], # yaw, roll, steer
                       [u4, u6, u7], # roll rate, rear wheel rate, steer rate
                       kd_eqs=kinematical,
                       q_dependent=[q5], # pitch angle
                       configuration_constraints=[holonomic],
                       u_dependent=[u3, u5, u8], # yaw rate, pitch rate, front wheel rate
                       velocity_constraints=nonholonomic)
fr, frstar = kane.kanes_equations(forces, bodies)

# Validation of non-linear equations

print('Loading numerical input parameters.')

from dtk import bicycle

bp = bicycle.benchmark_parameters()
mp = bicycle.benchmark_to_moore(bp)

# load the input values specified in table one of Basu-Mandal2007
basu_input = bicycle.basu_table_one_input()

# convert the values to my coordinates and speeds
moore_input = bicycle.basu_to_moore_input(basu_input, bp['rR'],
        bp['lam'])

constant_substitutions = {}
for k, v in mp.items():
    try:
        exec('constant_substitutions[{}] = v'.format(k))
    except NameError:
        pass
dynamic_substitutions = {}
for k, v in moore_input.items():
    try:
        exec('dynamic_substitutions[{}] = v'.format(k))
    except NameError:
        pass

substitutions = {T4: 0.0, T6: 0.0, T7: 0.0}
substitutions.update(constant_substitutions)
substitutions.update(dynamic_substitutions)

print('Substituting numerical parameters.')
mass_matrix = kane.mass_matrix_full.subs(substitutions)
forcing_vector = kane.forcing_full.subs(kane.kindiffdict()).subs(substitutions)

print("Solving for x'.")
xd = mass_matrix.inv() * forcing_vector

print("Generating output dictionary.")
# convert the outputs from my model to the Basu-Mandal coordinates
# TODO : raise an issue about not knowing which order the x vector is in with
# reference to M * x' = F
states = kane._q + kane._u
state_names = [str(state)[:-3] + 'p' for state in states]
moore_output = {k : v for k, v in zip(state_names, list(xd))}
u1 = -rr * u6 * sym.cos(q3)
u1p = u1.diff(t)
u2 = -rr * u6 * sym.sin(q3)
u2p = u2.diff(t)
moore_output['u1p'] = u1p.subs({u6d: moore_output['u6p']}).subs(kane.kindiffdict()).subs(substitutions)
moore_output['u2p'] = u2p.subs({u6d: moore_output['u6p']}).subs(kane.kindiffdict()).subs(substitutions)
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
        print('Failed: {} is supposed to be {} but is {}'.format(k, bv, mv))

# Try with lambdify

from pydy_code_gen.code import generate_ode_function
parameters = constant_substitutions.keys()
parameters.sort()
print('Generating a numeric right hand side function.')
rhs = generate_ode_function(kane, parameters, specified=[T4, T6, T7],
        generator="lambdify")
state_values = []
for state in kane._q + kane._u:
    state_values.append(substitutions[state])
parameter_values = [0.0, 0.0, 0.0] + [constant_substitutions[k] for k in sorted(constant_substitutions.keys())]
print('Evaluating the right hand side.')
xd = rhs(state_values, 0.0, parameter_values)

# yaw rate, and wheel contact diff equations
# kinetic and potential energy, and energy
# contact forces

#set_trace()

#def print_if(display, message):
    #if display == True:
        #print(message)

display = True

forces = {co: mc * g * N['3'],
          do: md * g * N['3'],
          eo: me * g * N['3'],
          fo: mf * g * N['3'],
          C: T4 * B['1'] - T6 * C['2'] - T7 * C['3'],
          D: T6 * C['2'],
          E: T7 * E['3']}

# TODO : probably need to sub in kind diffs somewhere in all this

# solve for the dependent u's: u3, u5, u8
# TODO : the following works but I need to make sure there are no remainder
# expressions, i.e. is nonholonomic[0] = c1 * u3 + ... + c6 * u8 + c7 where c7
# are terms that aren't linear in the u's
print("Solving for the dependent u's")
subs = {}
simple_holonomic = []
for equation in nonholonomic:
    new_equation = 0
    for speed in [u3, u4, u5, u6, u7, u8]:
        coefficient = equation.expand().coeff(speed)
        if coefficient != 0:
            z = sym.Dummy()
            subs[z] = coefficient
            new_equation += z * speed
    simple_holonomic.append(new_equation)

simple_dependent_u = sym.solve(simple_holonomic, u3, u5, u8, dict=True)[0]
dependent_u = {}
for key, equation in simple_dependent_u.items():
    dependent_u[key] = equation.subs(subs)
print("Dependent u's found.")

# differentiate the dependent u's
print("Differentiating the expressions for the dependent u's")
dependent_u_dot = {}
for k, v in dependent_u.items():
    dependent_u_dot[k.diff(t)] = v.diff(t)
print("U dots computed.")

# generalized active forces
print("Computing fr and frstar.")

fr = {}
frstar = {}
for u in [u4, u6, u7]:
    print('Computing the {} equations'.format(u))
    fr[u] = 0
    frstar[u] = 0
    for body in bodies:
        print('Adding the {} component'.format(body))

        try:
            R = forces[body.masscenter]
        except KeyError:
            R = 0

        try:
            T = forces[body.frame]
        except KeyError:
            T = 0

        vr = body.masscenter.vel(N).subs(kane.kindiffdict()).subs(dependent_u).diff(u, N)
        omega = body.frame.ang_vel_in(N).subs(kane.kindiffdict()).subs(dependent_u)
        wr = omega.diff(u, N)

        fr[u] += vr.dot(R) + wr.dot(T)

        a = body.masscenter.acc(N).subs(kane.kindiffdict()).subs(dependent_u).subs(dependent_u_dot)
        alpha = body.frame.ang_acc_in(N).subs(kane.kindiffdict()).subs(dependent_u).subs(dependent_u_dot)

        frstar[u] += (vr.dot(-body.mass * a) +
        wr.dot(-(alpha.dot(body.inertia[0]) +
            omega.cross(body.inertia[0]).dot(omega))))

u3p, u4p, u5p, u6p, u7p, u8p = sym.symbols('u3p u4p u5p u6p u7p u8p')
simple_u_dots = {u3d: u3p, u4d: u4p, u5d: u5p, u6d: u6p, u7d: u7p, u8d: u8p}

fr_equations = []
for u in [u4, u6, u7]:
    a = fr[u].subs(simple_u_dots).subs(substitutions)
    b = frstar[u].subs(simple_u_dots).subs(substitutions)
    fr_equations.append(a + b)

print("Done with fr and frstar.")
