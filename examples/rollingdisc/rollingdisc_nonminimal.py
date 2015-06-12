#!/usr/bin/env python

"""This script generates the equations of motion of a knife edged disc
rolling on a flat surface under the influence of gravity. This demonstration
is one of the simplest examples of a dynamic system which includes:

 - Configuration constraint
 - Velocity constraint
 - Auxiliary generalized forces to determine constraint forces

Keep in mind that the approach below is diabolical.  This problem can be
solved much more simply. The diabolical nature is only here to make sure
that the results from KanesMethod are correct.

"""

from sympy import symbols, Matrix, solve, Poly, sin, cos, zeros
from sympy.physics.mechanics import *

# Symbols for time and constant parameters
t, r, m, g, I, J = symbols('t r m g I J')
# Symbols for contact forces
Fx, Fy, Fz = symbols('Fx Fy Fz')

# Configuration variables and their time derivatives
# q[0] -- yaw
# q[1] -- lean
# q[2] -- spin
# q[3] -- dot(-r*B.z, A.z) -- distance from ground plane to disc center in A.z
#         direction
q = dynamicsymbols('q:4')
qd = [qi.diff(t) for qi in q]

# Generalized speeds and their time derivatives
# u[0] -- disc angular velocity component, disc fixed x direction
# u[1] -- disc angular velocity component, disc fixed y direction
# u[2] -- disc angular velocity component, disc fixed z direction
# u[3] -- disc velocity component, A.x direction
# u[4] -- disc velocity component, A.y direction
# u[5] -- disc velocity component, A.z direction
u = dynamicsymbols('u:6')
ud = [ui.diff(t) for ui in u]
ud_zero = {udi : 0 for udi in ud}

# Auxiliary generalized speeds
# ua[0] -- contact point auxiliary generalized speed, A.x direction
# ua[1] -- contact point auxiliary generalized speed, A.y direction
# ua[2] -- contact point auxiliary generalized speed, A.z direction
ua = dynamicsymbols('ua:3')
ua_zero = {uai : 0 for uai in ua}

# Reference frames
N = ReferenceFrame('N')
A = N.orientnew('A', 'Axis', [q[0], N.z])   # Yaw intermediate frame
B = A.orientnew('B', 'Axis', [q[1], A.x])   # Lean intermediate frame
C = B.orientnew('C', 'Axis', [q[2], B.y])   # Disc fixed frame

# Angular velocity and angular acceleration of disc fixed frame
C.set_ang_vel(N, u[0]*B.x + u[1]*B.y + u[2]*B.z)
C.set_ang_acc(N, C.ang_vel_in(N).diff(t, B)
               + cross(B.ang_vel_in(N), C.ang_vel_in(N)))

# Velocity and acceleration of points
P = Point('P')                   # Disc-ground contact point
P.set_vel(N, ua[0]*A.x + ua[1]*A.y + ua[2]*A.z)
O = P.locatenew('O', q[3]*A.z + r*sin(q[1])*A.y)     # Center of disc
O.set_vel(N, u[3]*A.x + u[4]*A.y + u[5]*A.z)
O.set_acc(N, O.vel(N).diff(t, A) + cross(A.ang_vel_in(N), O.vel(N)))

# Kinematic differential equations
w_c_n_qd = qd[0]*A.z + qd[1]*B.x + qd[2]*B.y
v_o_n_qd = O.pos_from(P).diff(t, A) + cross(A.ang_vel_in(N), O.pos_from(P))
kindiffs = Matrix([dot(w_c_n_qd - C.ang_vel_in(N), uv) for uv in B] +
                  [dot(v_o_n_qd - O.vel(N), A.z)])
qd_kd = solve(kindiffs, qd)     # solve for dq/dt's in terms of u's
print("Kinematic differential equations")
mprint(qd_kd)

# Values of generalized speeds during a steady turn
steady_conditions = solve(kindiffs.subs({qd[1] : 0, qd[3] : 0}), u)
steady_conditions.update({qd[1] : 0, qd[3] : 0})
print("Steady turning conditions")
mprint(steady_conditions)

# Partial angular velocities and velocities
partial_w_C = [C.ang_vel_in(N).diff(ui, N) for ui in u + ua]
partial_v_O = [O.vel(N).diff(ui, N) for ui in u + ua]
partial_v_P = [P.vel(N).diff(ui, N) for ui in u + ua]

print("Partial angular velocities")
mprint(partial_w_C)
print("Partial velocities of O")
mprint(partial_v_O)
print("Partial velocities of P")
mprint(partial_v_P)

# Configuration constraint
print("Configuration constraint")
f_c = Matrix([dot(-r*B.z, A.z) - q[3]])
# Velocity constraints
print("Velocity constraints")
f_v = Matrix([dot(O.vel(N) - (P.vel(N) + cross(C.ang_vel_in(N),
    O.pos_from(P))), ai).expand() for ai in A])
mprint(f_v)
# Acceleration constraints
v_o_n = cross(C.ang_vel_in(N), O.pos_from(P))
a_o_n = v_o_n.diff(t, A) + cross(A.ang_vel_in(N), v_o_n)
f_a = Matrix([dot(O.acc(N) - a_o_n, ai) for ai in A])
print("Acceleration constraints")
mprint(f_a)

# Constraint coefficient matrix:  M_v * [u; ua] = 0
M_v = zeros(3, 9)
for i in range(3):
    for j, ui in enumerate(u + ua):
        M_v[i, j] = f_v[i].diff(ui)

# Dependent speeds in terms of independent and auxiliary speeds
# After this point, we have married ourselves this choice of independent speeds
M_v_i = M_v[:, :3]       # take u[0], u[1], u[2] as independent
M_v_d = M_v[:, 3:6]      # take u[3], u[4], u[5] as dependent
M_v_aux = M_v[:, 6:]     # last colums are the auxiliary
M_v_i_aux = M_v_i.row_join(M_v_aux)
A_rs = - M_v_d.inv() * M_v_i_aux
# Constraint matrix:  u_dependent = A_rs * [u_i; u_aux]
mprint(A_rs)

# Dependent speeds under the condition u_aux == 0
u_dep = A_rs[:, :3] * Matrix(u[:3])
u_dep_dict = {udi : u_depi[0] for udi, u_depi in zip(u[3:], u_dep.tolist())}
mprint(u_dep_dict)

# Active forces
F_O = m*g*A.z
F_P = Fx * A.x + Fy * A.y + Fz * A.z
# Generalized active forces (unconstrained)
Fr_u = Matrix([dot(F_O, pv_o) + dot(F_P, pv_p) for pv_o, pv_p in
        zip(partial_v_O, partial_v_P)])


# Inertia force
R_star_O = -m*O.acc(N)

# Inertia torque
I_C_O = inertia(B, I, J, I)
T_star_C = -(dot(I_C_O, C.ang_acc_in(N)) \
             + cross(C.ang_vel_in(N), dot(I_C_O, C.ang_vel_in(N))))

# Generalized inertia forces (unconstrained)
Fr_star_u = Matrix([dot(R_star_O, pv) + dot(T_star_C, pav) for pv, pav in
                    zip(partial_v_O, partial_w_C)])

# Form nonholonomic Fr and nonholonomic Fr_star
# See equations 4.4.3 and  4.11.4 of Kane & Levinson
Fr_c = Fr_u[:3, :].col_join(Fr_u[6:, :]) + A_rs.T * Fr_u[3:6, :]
Fr_star_c = Fr_star_u[:3, :].col_join(Fr_star_u[6:, :])\
            + A_rs.T * Fr_star_u[3:6, :]
Fr_star_steady = Fr_star_c.subs(ud_zero).subs(u_dep_dict)\
        .subs(steady_conditions).subs({q[3]: -r*cos(q[1])}).expand()

mprint(Fr_c)
mprint(Fr_star_steady)

# First dynamic equation, under steady conditions is 2nd order polynomial in
# dq0/dt.
steady_turning_dynamic_equation = Fr_c[0] + Fr_star_steady[0]
# Equilibrium is posible when the solution to this quadratic is real, i.e.,
# when the discriminant in the quadratic is non-negative
p = Poly(steady_turning_dynamic_equation, qd[0])
a, b, c = p.coeffs()
discriminant = b*b - 4*a*c      # Must be non-negative for equilibrium
# in case of thin disc inertia assumptions
#mprint((discriminant / (r**3 * m**2)).expand())


# ADD ALL CODE DIRECTLY BELOW HERE, do not change above!
# Think there should be at 12 assertion tests:
# 1) Fr[i] == fr from KanesMethod  i = 0, ..., 5
# 2) Fr_star[i] == frstar from KanesMethod i = 0, ..., 5
# if 2) is slow, try comparing this instead:
# 2a) Fr_star_steady[i] == frstar from KanesMethod, evaluated at steady turning
# conditions.
# This should be something like frstar.subs(ud_zero).subs(steady_conditions)
