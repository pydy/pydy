"""Rolling disc with a trailing/leading skate attached by a massless arm."""
from sympy import symbols, Matrix, solve, Poly, trigsimp, zeros
from sympy.physics.mechanics import (cross, dot, dynamicsymbols, inertia,
     ReferenceFrame, Point)
# Symbols for time and constant parameters
t, r, l, m, g, tau, F_friction = symbols('t r l m g tau F_friction')
# Configuration variables and their time derivatives
# q[0] -- yaw
# q[1] -- lean
# q[2] -- spin
# q[3] -- arm angle relative to line from contact point to disc center
q = dynamicsymbols('q:4')
qd = [qi.diff(t) for qi in q]
# Generalized speeds and their time derivatives
# u[0] -- disc angular velocity component, roll frame x direction
# u[1] -- disc angular velocity component, roll frame y direction
# u[2] -- disc angular velocity component, roll frame z direction
# u[3] -- arm angular velocity relative to roll frame
u = dynamicsymbols('u:4')
ud = [ui.diff(t) for ui in u]
ud_zero = {udi : 0 for udi in ud}
## Orientation, angular velocity, and, angular acceleration of frames
# Orientation of all frames
N = ReferenceFrame('N')			    # Inertial frame, N.z down
A = N.orientnew('A', 'Axis', [q[0], N.z])   # Yaw frame
B = A.orientnew('B', 'Axis', [q[1], A.x])   # Roll frame
C = B.orientnew('C', 'Axis', [q[2], B.y])   # Disc fixed frame
D = B.orientnew('D', 'Axis', [q[3], B.y])   # Arm fixed frame
# Angular velocity of disc fixed frame and arm fixed frame
C.set_ang_vel(N, u[0] * B.x + u[1] * B.y + u[2] * B.z)
D.set_ang_vel(N, C.ang_vel_in(N) + u[3] * B.y)
omega_C_N_qd = qd[0] * A.z + qd[1] * B.x + qd[2] * B.y
omega_D_N_qd = qd[0] * A.z + qd[1] * B.x + (qd[2] + qd[3]) * B.y
# Form the kinematic ODE's which relate dq/dt with u
kinematic_odes = ([dot(omega_C_N_qd - C.ang_vel_in(N), uv) for uv in B]
		  + [dot(omega_D_N_qd - D.ang_vel_in(N), B.y)])
qdots = solve(kinematic_odes, qd)
for qdi in qd:
    print("{0} = {1}".format(qdi, qdots[qdi]))
# Eliminate dq/dt from angular velocity expressions
A.set_ang_vel(N, A.ang_vel_in(N).subs(qdots))
B.set_ang_vel(N, B.ang_vel_in(N).subs(qdots))
D.set_ang_vel(N, D.ang_vel_in(N).subs(qdots))
print("omega_A_N = {0}".format(A.ang_vel_in(N)))
print("omega_B_N = {0}".format(B.ang_vel_in(N)))
print("omega_C_N = {0}".format(C.ang_vel_in(N)))
print("omega_D_N = {0}".format(D.ang_vel_in(N)))
# Angular acceleration of disc and arm
C.set_ang_acc(N, (C.ang_vel_in(N).diff(t, B)
                 + cross(B.ang_vel_in(N), C.ang_vel_in(N))).subs(qdots))
D.set_ang_acc(N, (D.ang_vel_in(N).diff(t, B)
	         + cross(B.ang_vel_in(N), D.ang_vel_in(N))).subs(qdots))
print("alpha_C_N = {0}".format(C.ang_acc_in(N)))
print("alpha_D_N = {0}".format(D.ang_acc_in(N)))
# Position, velocity, and, acceleration of points
O = Point('O')                   # Disc-ground contact point
P = O.locatenew('P', -r*B.z)     # Center of disc
Q = P.locatenew('Q', l*D.z)      # Arm mass center, arm-ground contact point
# Velocity and acceleration of disc center
P.set_vel(N, cross(C.ang_vel_in(N), P.pos_from(O)))
P.set_acc(N, (P.vel(N).diff(t, B)
              + cross(B.ang_vel_in(N), P.vel(N))).subs(qdots))
print("v_P_N = {0}".format(P.vel(N)))
print("a_P_N = {0}".format(P.acc(N)))
# Velocity and acceleration of arm-ground contact point
Q.v2pt_theory(P, N, D)
Q.a2pt_theory(P, N, D)
print("v_Q_N = {0}".format(Q.vel(N)))
print("a_Q_N = {0}".format(Q.acc(N)))
# Configuration constraint
f_c = dot(Q.pos_from(O), A.z)
print("f_c = {0}".format(f_c))
# Velocity constraints
# No lateral slip of arm and no velocity perpendicular to ground plane
#f_v = [dot(Q.vel(N), A.y).expand(), f_c.diff(t).subs(qdots).expand()]
f_v = [dot(Q.vel(N), A.y).expand(), dot(Q.vel(N), A.z)]
print("f_v[0] = {0}".format(f_v[0]))
print("f_v[1] = {0}".format(f_v[1]))
f_v_du = zeros((2, 4))
for i in range(2):
    for j in range(4):
        f_v_du[i, j] = f_v[i].diff(u[j])
print("Velocity constraint matrix: ")
print(f_v_du)
# Matrix mapping inpdependent speeds to dependent speeds
A_rs = trigsimp((-f_v_du[:, 2:].inv() * f_v_du[:, :2]).expand())
print("A_rs = \n{0}".format(A_rs))
# Steady turning conditions generalized speeds during a steady turn
qd_steady = {qd[1] : 0, qd[3] : 0}
kinematic_odes_steady = [eq.subs(qd_steady) for eq in kinematic_odes]
u_steady = solve(kinematic_odes_steady, u)
for ui in u:
    print("{0}_steady = {1}".format(ui, u_steady[ui]))
# Partial angular velocities and velocities
partial_w_C_N = [C.ang_vel_in(N).diff(ui, N) for ui in u]
partial_w_D_N = [D.ang_vel_in(N).diff(ui, N) for ui in u]
partial_v_P_N = [P.vel(N).diff(ui, N) for ui in u]
partial_v_Q_N = [Q.vel(N).diff(ui, N) for ui in u]
print("Partial angular velocities of C in N: {0}".format(partial_w_C_N))
print("Partial angular velocities of D in N: {0}".format(partial_w_D_N))
print("Partial velocities of P in N: {0}".format(partial_v_P_N))
print("Partial velocities of Q in N: {0}".format(partial_v_Q_N))
# Active torques
T_C = tau * B.y
T_D = -tau * B.y
# Active Forces
F_P = m * g * A.z
F_Q = F_friction * A.x
# Generalized active forces
F = Matrix([dot(T_C, pav_c)
          + dot(T_D, pav_d)
          + dot(F_P, pv_p)
          + dot(F_Q, pv_q)
          for pav_c, pav_d, pv_p, pv_q in
          zip(partial_w_C_N, partial_w_D_N, partial_v_P_N, partial_v_Q_N)])
print("Generalized active forces:")
for i, F_i in enumerate(F):
    print("F[{0}] = {1}".format(i, F_i))
F_constrained = F[:2, 0] + A_rs.transpose() * F[2:, 0]
print("Constrained generalized active forces:")
for i, F_i in enumerate(F_constrained):
    print("F[{0}] = {1}".format(i, F_i))
# Inertia torque
I_C_P = inertia(B, m*r*r/4, m*r*r/2, m*r*r/4)
T_star_C_N = -(dot(I_C_P, C.ang_acc_in(N))
               + cross(C.ang_vel_in(N), dot(I_C_P, C.ang_vel_in(N))))
# Inertia forces
R_star_P_N = -m * P.acc(N)
# Generalized inertia forces
F_star = Matrix([(dot(T_star_C_N, pav_c)
                + dot(R_star_P_N, pv_p)).simplify()
                for pav_c, pv_p, pv_q in
                zip(partial_w_C_N, partial_v_P_N, partial_v_Q_N)])
print("Generalized inertia forces:")
for i, F_star_i in enumerate(F_star):
    print("F_star[{0}] = {1}".format(i, F_star_i))
F_star_constrained = F_star[:2, 0] + A_rs.transpose() * F_star[2:, 0]
print("Constrained generalized inertia forces:")
for i, F_star_i in enumerate(F_star_constrained):
    print("F_star[{0}] = {1}".format(i, F_star_i))
F_star_steady = [F_star_i.subs(ud_zero).subs(u_steady).expand().simplify()
                 for F_star_i in F_star_constrained]
print("Constrained generalized inertia forces (steady):")
for i, F_star_steady_i in enumerate(F_star_steady):
    print("F_star_steady[{0}] = {1}".format(i, F_star_steady_i))
# Under steady turning equations, the first dynamic equation is a quadratic in
# dq0/dt.  Equilibrium is possible when the solution to this quadratic is real,
# i.e., when the discriminant is non-negative
steady_turning_equation = ((F_constrained[0] + F_star_steady[0])/m/r/r).expand()
p = Poly(steady_turning_equation, qd[0])
a, b, c = p.coeffs()
print("a = {0}".format(a))
print("b = {0}".format(b))
print("c = {0}".format(c))
discriminant = b*b - 4*a*c      # Must be non-negative for equilibrium
print("discriminant = {0}".format(discriminant))
q2_dot = solve(discriminant, qd[2])
# Boundary of steady turning equilibria; to one side of this line solutions are
# impossible.
print(q2_dot[0].simplify())
print(q2_dot[1].simplify())
