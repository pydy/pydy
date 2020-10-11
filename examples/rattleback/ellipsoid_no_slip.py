from sympy import (symbols, sqrt, ccode, acos, Symbol, sin,
    cos, tan, cse, numbered_symbols)
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Vector,
    Point, inertia, dot, cross)

Vector.simp = False         # Prevent the use of trigsimp and simplify
t, g = symbols('t g')       # Time and gravitational constant
a, b, c = symbols('a b c')  # semi diameters of ellipsoid
d, e, f = symbols('d e f')  # mass center location parameters
s = symbols('s')            # coefficient of viscous friction

# Mass and Inertia scalars
m, Ixx, Iyy, Izz, Ixy, Iyz, Ixz = symbols('m Ixx Iyy Izz Ixy Iyz Ixz')

q = dynamicsymbols('q:5')             # Generalized coordinates
qd = [qi.diff(t) for qi in q]         # Generalized coordinate time derivatives
u = dynamicsymbols('u:3')             # Generalized speeds
ud = [ui.diff(t) for ui in u]         # Generalized speeds time derivatives
ua = dynamicsymbols('ua:3')           # Auxiliary generalized speeds
CF = dynamicsymbols('Rx Ry Rz')       # Contact forces
r = dynamicsymbols('r:3')             # Coordinates, in R frame, from O to P
rd = [ri.diff(t) for ri in r]         # Time derivatives of xi

N = ReferenceFrame('N')                   # Inertial Reference Frame
Y = N.orientnew('Y', 'Axis', [q[0], N.z]) # Yaw Frame
L = Y.orientnew('L', 'Axis', [q[1], Y.x]) # Lean Frame
R = L.orientnew('R', 'Axis', [q[2], L.y]) # Rattleback body fixed frame

I = inertia(R, Ixx, Iyy, Izz, Ixy, Iyz, Ixz)    # Inertia dyadic

# Angular velocity using u's as body fixed measure numbers of angular velocity
R.set_ang_vel(N, u[0]*R.x + u[1]*R.y + u[2]*R.z)

# Rattleback ground contact point
P = Point('P')
P.set_vel(N, ua[0]*Y.x + ua[1]*Y.y + ua[2]*Y.z)

# Rattleback ellipsoid center location, see:
# "Realistic mathematical modeling of the rattleback", Kane, Thomas R. and
# David A. Levinson, 1982, International Journal of Non-Linear Mechanics
mu = [dot(rk, Y.z) for rk in R]
eps = sqrt((a*mu[0])**2 + (b*mu[1])**2 + (c*mu[2])**2)
O = P.locatenew('O', -a*a*mu[0]/eps*R.x
                     -b*b*mu[1]/eps*R.y
                     -c*c*mu[2]/eps*R.z)
O.set_vel(N, P.vel(N) + cross(R.ang_vel_in(N), O.pos_from(P)))

# Mass center position and velocity
RO = O.locatenew('RO', d*R.x + e*R.y + f*R.z)
RO.set_vel(N, P.vel(N) + cross(R.ang_vel_in(N), RO.pos_from(P)))

qd_rhs = [(u[2]*cos(q[2]) - u[0]*sin(q[2]))/cos(q[1]),
         u[0]*cos(q[2]) + u[2]*sin(q[2]),
         u[1] + tan(q[1])*(u[0]*sin(q[2]) - u[2]*cos(q[2])),
         dot(P.pos_from(O).diff(t, R), N.x),
         dot(P.pos_from(O).diff(t, R), N.y)]

# Partial angular velocities and partial velocities
partial_w = [R.ang_vel_in(N).diff(ui, N) for ui in u + ua]
partial_v_P = [P.vel(N).diff(ui, N) for ui in u + ua]
partial_v_RO = [RO.vel(N).diff(ui, N) for ui in u + ua]

# Set auxiliary generalized speeds to zero in velocity vectors
P.set_vel(N, P.vel(N).subs({ua[0]:0, ua[1]:0, ua[2]:0}))
O.set_vel(N, O.vel(N).subs({ua[0]:0, ua[1]:0, ua[2]:0}))
RO.set_vel(N, RO.vel(N).subs({ua[0]:0, ua[1]:0, ua[2]:0}))

# Angular acceleration
R.set_ang_acc(N, ud[0]*R.x + ud[1]*R.y + ud[2]*R.z)

# Acceleration of mass center
RO.set_acc(N, RO.vel(N).diff(t, R) + cross(R.ang_vel_in(N), RO.vel(N)))

# Forces and Torques
F_P = sum([cf*uv for cf, uv in zip(CF, Y)], Vector(0))
F_RO = m*g*Y.z
T_R = -s*R.ang_vel_in(N)

# Generalized Active forces
gaf_P = [dot(F_P, pv) for pv in partial_v_P]
gaf_RO = [dot(F_RO, pv) for pv in partial_v_RO]
gaf_R = [dot(T_R, pv) for pv in partial_w]

# Generalized Inertia forces
# First, compute R^* and T^* for the rigid body
R_star = - m*RO.acc(N)
T_star = - dot(R.ang_acc_in(N), I)\
         - cross(R.ang_vel_in(N), dot(I, R.ang_vel_in(N)))

# Isolate the parts that involve only time derivatives of u's
R_star_udot = sum([R_star.diff(udi, N)*udi for udi in ud], Vector(0))
T_star_udot = sum([T_star.diff(udi, N)*udi for udi in ud], Vector(0))
for ui in u:
  assert(R_star_udot.diff(ui, N) == 0)
  assert(T_star_udot.diff(ui, N) == 0)

# Isolate the parts that involve no time derivatives of u's
R_star_other = R_star.subs({ud[0]:0, ud[1]:0, ud[2]:0})
T_star_other = T_star.subs({ud[0]:0, ud[1]:0, ud[2]:0})
for udi in ud:
  assert(R_star_other.diff(udi, N) == 0)
  assert(T_star_other.diff(udi, N) == 0)

gif_udot = []     # Generalized inertia forces with udots
gif_other = []    # Generalized inertia forces without udots
for i in range(len(u + ua)):
  gif_udot.append(dot(partial_w[i], T_star_udot)
                + dot(partial_v_RO[i], R_star_udot))
  gif_other.append(dot(partial_w[i], T_star_other)
                 + dot(partial_v_RO[i], R_star_other))

# The first three equations of Fr + Fr^* = 0 are the dynamic equations
# associated with the three independent generalized speeds, u0, u1, u2.  These
# equations ultimately need to be solved for the time derivatives of the u's,
# so with this in mind, we rearrange them as:
# M_dyn(q) * du/dt = f_dyn(q, u)
f_dyn = [0]*3
M_dyn = [0]*9
for i in range(3):
  f_dyn[i] = - (gaf_P[i] + gaf_RO[i] + gaf_R[i] + gif_other[i])
  for j in range(3):
    M_dyn[3*i + j] = gif_udot[i].diff(ud[j])

# The last three equations of Fr + Fr^* = 0 are the auxiliary dynamic equations
# associated with the three auxiliary generalized speeds.  These equations
# ultimately need to be solved for the constraint forces.  With this in mind we
# rearrange them as:
# CF = f_cf(q, u, ud)
f_cf = [0]*3
for i in range(3):
  f_cf[i] = - (gaf_RO[i + 3] + gaf_R[i + 3] + gif_udot[i + 3] + gif_other[i + 3])
  assert(gaf_P[i + 3] == CF[i])

# Kinetic and potential energy
ke = (m*dot(RO.vel(N), RO.vel(N)) + dot(R.ang_vel_in(N), dot(I, R.ang_vel_in(N))))/2.0
pe = -m*g*dot(RO.pos_from(P), Y.z)

# Delta -- angle between Y.z and R.z
delta = acos(dot(Y.z, R.z))

# Jacobian matrix
J = [0]*25    # only consider the 5 states (ignore x, y, yaw ode's)
for i, de_rhs in enumerate(qd_rhs[1:3] + f_dyn):
  for j, xi in enumerate(q[1:3] + u):
    J[5*i + j] = de_rhs.diff(xi)
    for qdk, qdk_rhs in zip(qd, qd_rhs):
      J[5*i + j] += de_rhs.diff(qdk)*qdk_rhs.diff(xi)

# Build lists of grouped equations to do CSE on
exp_ode = qd_rhs + M_dyn + f_dyn
exp_output = f_cf + [ke, pe, ke + pe, delta]
exp_jac = J + M_dyn

# Subsitution dictionary to replace dynamic symbols with regular symbols
subs_dict = {q[0]: Symbol('q0'), q[1]: Symbol('q1'), q[2]: Symbol('q2'),
             qd[0]: Symbol('qd0'), qd[1]: Symbol('qd1'), qd[2]: Symbol('qd2'),
             u[0]: Symbol('u0'), u[1]: Symbol('u1'), u[2]: Symbol('u2'),
             ud[0]: Symbol('ud0'), ud[1]: Symbol('ud1'), ud[2]: Symbol('ud2')}

for i in range(len(exp_ode)):
  exp_ode[i] = exp_ode[i].subs(subs_dict)

for i in range(len(exp_output)):
  exp_output[i] = exp_output[i].subs(subs_dict)

for i in range(len(exp_jac)):
  exp_jac[i] = exp_jac[i].subs(subs_dict)

# CSE on all quantities needed for numerical integration of ordinary
# differential equations:  qd_rhs (5), M_dyn (9), f_dyn (3)
z, exp_ode_red = cse(exp_ode, numbered_symbols('z'))

output_code = "  // Intermediate variables for ODE function\n"
for zi_lhs, zi_rhs in z:
  output_code += "  {0} = {1};\n".format(zi_lhs, ccode(zi_rhs))

output_code += "\n  // Kinematic differential equations\n"
for i in range(5):
  output_code += "  dxdt[{0}] = {1};\n".format(i, ccode(exp_ode_red[i]))

output_code += "\n  // Mass matrix\n"
for i in range(3):
  for j in range(3):
    output_code += "  M_dyn({0}, {1}) = {2};\n".format(i, j,
                                ccode(exp_ode_red[5 + 3*i + j]))

output_code += "\n  // Right hand side of dynamic equations\n"
for i in range(3):
  output_code += "  f_dyn({0}) = {1};\n".format(i, ccode(exp_ode_red[14 + i]))

# CSE on all output quantities:  CF (3), ke, pe, te, delta
output_code += "\n  // Output quantites (evaluated at each output time-step)\n"
z, exp_output_red = cse(exp_output, numbered_symbols('z'))
for zi_lhs, zi_rhs in z:
  output_code += "  {0} = {1};\n".format(zi_lhs, ccode(zi_rhs))

output_code += "\n  // Contact forces\n"
for i in range(3):
  output_code += "  sd->CF[{0}] = {1};\n".format(i, ccode(exp_output_red[i]))

output_code += "\n  // Mechanical energy\n"
for i, name in enumerate(["ke", "pe", "te"]):
  output_code += "  sd->{0} = {1};\n".format(name,
                                            ccode(exp_output_red[i + 3]))
output_code += "  // Tilt of Rattleback with respect to vertical\n"
output_code += "  sd->delta = {0};\n".format(ccode(exp_output_red[-1]))

# CSE on all quantities needed for Jacobian matrix:  M_dyn (3), J (64)
output_code += "\n  // Intermediate quantities needed for Jacobian matrix\n"
z, exp_jac_red = cse(exp_jac, numbered_symbols('z'))
for zi_lhs, zi_rhs in z:
  output_code += "  " + str(zi_lhs) + " = " + ccode(zi_rhs) + ";\n"

output_code += "\n  // Entries of Jacobian matrix\n"
for i in range(5):
  for j in range(5):
    output_code += "  J({0}, {1}) = {2};\n".format(i, j,
        ccode(exp_jac_red[5*i + j]))
output_code += "\n  // Entries of Mass matrix\n"
for i in range(3):
  for j in range(3):
    output_code += "  M_dyn({0}, {1}) = {2};\n".format(i, j,
        ccode(exp_jac_red[25 + 3*i + j]))

# Perform text substitutions to change symbols used for state variables and
# their derivatives (qi, ui, qdi, udi) to the names used by the ode integrator.
import re
output_code = re.sub(r"z(\d+)", r"z[\1]", output_code)
output_code = re.sub(r"q([01234])", r"x[\1]", output_code)
output_code = re.sub(r"qd([01234])", r"dxdt[\1]", output_code)
output_code = re.sub(r"u([012])", r"x[\1 + 5]", output_code)
output_code = re.sub(r"ud([012])", r"dxdt[\1 + 5]", output_code)

with open("ellipsoid_no_slip.txt", 'w') as f:
    f.write(output_code)
