from sympy import symbols, Matrix, solve, Poly
from sympy.physics.mechanics import *

# Symbols for time and constant parameters
t, r, m, g, I, J = symbols('t r m g I J')
# Symbols for contact forces
Fx, Fy, Fz = symbols('Fx Fy Fz')

# Configuration variables and their time derivatives
# q[0] -- yaw
# q[1] -- lean
# q[2] -- spin
q = dynamicsymbols('q:3')
qd = [qi.diff(t) for qi in q]

# Generalized speeds and their time derivatives
# u[0] -- disc angular velocity component, disc fixed x direction
# u[1] -- disc angular velocity component, disc fixed y direction
# u[2] -- disc angular velocity component, disc fixed z direction
u = dynamicsymbols('u:3')
ud = [ui.diff(t) for ui in u]
ud_zero = {udi : 0 for udi in ud}     #

# Auxiliary generalized speeds
# ua[0] -- contact point auxiliary generalized speed, x direction
# ua[1] -- contact point auxiliary generalized speed, y direction
# ua[2] -- contact point auxiliary generalized speed, z direction
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
O = P.locatenew('O', -r*B.z)     # Center of disc
P.set_vel(N, ua[0]*A.x + ua[1]*A.y + ua[2]*A.z)
O.v2pt_theory(P, N, C)
O.set_acc(N, O.vel(N).subs(ua_zero).diff(t, B)
           + cross(B.ang_vel_in(N), O.vel(N).subs(ua_zero)))

# Kinematic differential equations
w_c_n_qd = qd[0]*A.z + qd[1]*B.x + qd[2]*B.y
kindiffs = Matrix([dot(w_c_n_qd - C.ang_vel_in(N), uv) for uv in B])
qd_kd = solve(kindiffs, qd)     # solve for dq/dt's in terms of u's
mprint(kindiffs)

# Values of generalized speeds during a steady turn
steady_conditions = solve(kindiffs.subs({qd[1] : 0}), u)
steady_conditions.update({qd[1] : 0})
print(steady_conditions)

# Partial angular velocities and velocities
partial_w_C = [C.ang_vel_in(N).diff(ui, N) for ui in u + ua]
partial_v_O = [O.vel(N).diff(ui, N) for ui in u + ua]
partial_v_P = [P.vel(N).diff(ui, N) for ui in u + ua]

print(partial_w_C)
print(partial_v_O)
print(partial_v_P)

# Active forces
F_O = m*g*A.z
F_P = Fx * A.x + Fy * A.y + Fz * A.z
# Generalized active forces
Fr = [dot(F_O, pv_o) + dot(F_P, pv_p) for pv_o, pv_p in
        zip(partial_v_O, partial_v_P)]

# Inertia force
R_star_O = -m*O.acc(N)

# Inertia torque
I_C_O = inertia(B, I, J, I)
T_star_C = -(dot(I_C_O, C.ang_acc_in(N)) \
             + cross(C.ang_vel_in(N), dot(I_C_O, C.ang_vel_in(N))))

# Generalized inertia forces
Fr_star = [dot(R_star_O, pv) + dot(T_star_C, pav) for pv, pav in
           zip(partial_v_O, partial_w_C)]


Fr_star_steady = [Fr_star_i.subs(ud_zero).subs(steady_conditions).expand()
                  for Fr_star_i in Fr_star]

mprint(Fr)
mprint(Fr_star_steady)

# First dynamic equation, under steady conditions is 2nd order polynomial in
# dq0/dt.
steady_turning_dynamic_equation = Fr[0] + Fr_star_steady[0]
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
