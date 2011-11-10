from sympy import (symbols, sqrt, zeros, ccode, acos, Symbol, sin,
    cos, tan)
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Vector,
    Point, inertia, dot, cross)

Vector.simp = False         # Prevent the use of trigsimp and simplify
t, g = symbols('t g')        # Time and gravitational constant
a, b, c = symbols('a b c')  # semi diameters of ellipsoid
d, e, f = symbols('d e f')  # mass center location parameters

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
#mu = [dot(rk, Y.z) for rk in R]
#eps = sqrt((a*mu[0])**2 + (b*mu[1])**2 + (c*mu[2])**2)
#O = P.locatenew('O', -a*a*mu[0]/eps*R.x
#                     -b*b*mu[1]/eps*R.y
#                     -c*c*mu[2]/eps*R.z)
O = P.locatenew('O', sum([-ri*uv for ri, uv in zip(r, R)]))
O.set_vel(N, P.vel(N) + cross(R.ang_vel_in(N), O.pos_from(P)))

# Mass center position, velocity, acceleration
RO = O.locatenew('RO', d*R.x + e*R.y + f*R.z)
RO.set_vel(N, P.vel(N) + cross(R.ang_vel_in(N), RO.pos_from(P)))

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

# Forces
F_P = sum([cf*uv for cf, uv in zip(CF, Y)])
F_RO = m*g*Y.z

# Generalized Active forces
gaf_P = [dot(F_P, pv) for pv in partial_v_P]
gaf_RO = [dot(F_RO, pv) for pv in partial_v_RO]

# Generalized Inertia forces
# First, compute R^* and T^* for the rigid body
R_star = - m*RO.acc(N)
T_star = - dot(R.ang_acc_in(N), I)\
         - cross(R.ang_vel_in(N), dot(I, R.ang_vel_in(N)))

# Isolate the parts that involve only time derivatives of u's
R_star_udot = sum([R_star.diff(udi, N)*udi for udi in ud])
T_star_udot = sum([T_star.diff(udi, N)*udi for udi in ud])
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
f_dyn = zeros(3, 1)
M_dyn = zeros(3, 3)
for i in range(3):
  f_dyn[i, 0] = - (gaf_P[i] + gaf_RO[i] + gif_other[i])
  for j in range(3):
    M_dyn[i, j] = gif_udot[i].diff(ud[j])
udots = M_dyn.LUsolve(f_dyn).T.tolist()[0]

# The last three equations of Fr + Fr^* = 0 are the auxiliary dynamic equations
# associated with the three auxiliary generalized speeds.  These equations
# ultimately need to be solved for the constraint forces.  With this in mind we
# rearrange them as:
# CF = f_cf(q, u, ud)
f_cf = zeros(3, 1)
for i in range(3):
  f_cf[i] = - (gaf_RO[i + 3] + gif_udot[i + 3] + gif_other[i + 3])
  assert(gaf_P[i + 3] == CF[i])

# Kinematic differential equations
# Angular coordinates:
qdots = [(u[2]*cos(q[2]) - u[0]*sin(q[2]))/cos(q[1]),
         u[0]*cos(q[2]) + u[2]*sin(q[2]),
         u[1] + tan(q[1])*(u[0]*sin(q[2]) - u[2]*cos(q[2])),
         dot(P.pos_from(O).diff(t, R), N.x),
         dot(P.pos_from(O).diff(t, R), N.y)]

# Jacobian matrix
# TODO: add partials of the position vector measure numbers as well as the
# partial derivatives of the position vector measure number time derivatives
J = zeros(8, 8)
for de_rhs in qdots + udots:
  for xi in q + u:
    J[i, j] = de_rhs.diff(xi)

# Kinetic and potential energy
ke = (m*dot(RO.vel(N), RO.vel(N)) + dot(R.ang_vel_in(N), dot(I, R.ang_vel_in(N))))/2.0
pe = -m*g*dot(RO.pos_from(P), Y.z)

# Delta -- angle between Y.z and R.z
delta = acos(dot(Y.z, R.z))

subs_dict = {qd[0]: Symbol('qd0'), qd[1]: Symbol('qd1'), qd[2]: Symbol('qd2'),
             ud[0]: Symbol('ud0'), ud[1]: Symbol('ud1'), ud[2]: Symbol('ud2'),
             r[0]: Symbol('r0'), r[1]: Symbol('r1'), r[2]: Symbol('r2'),
             rd[0]: Symbol('rd0'), rd[1]: Symbol('rd1'), rd[2]: Symbol('rd2')}

output_qdots = ""
for i, qdi in enumerate(qd):
  output_qdots += "  dxdt[" + str(i) + "] = " + ccode(qdots[i].subs(subs_dict)) + ";\n"

output_udots = ""
for i, udi in enumerate(ud):
  output_udots += "  dxdt[" + str(i+5) + "] = " + ccode(udots[i].subs(subs_dict)) + ";\n"

# We output a few extra quantities at each time step:
# -- contact forces
# -- kinetic, potential, and total energy
# -- delta
output_extra = ""
for i, cfi in enumerate(CF):
  output_extra += "  CF[" + str(i) + "] = " + ccode(f_cf[i].subs(subs_dict)) + ";\n"
output_extra += "  ke = " + ccode(ke) + ";\n"
output_extra += "  pe = " + ccode(pe) + ";\n"
output_extra += "  te = ke + pe;\n"
output_extra += "  delta = " + ccode(delta) + ";\n"

output = (output_qdots + "\n\n" + output_udots + "\n\n" + output_extra).replace("(t)", "")

state_rep = [("q{0}".format(i), "x[{0}]".format(i)) for i in range(5)] +\
            [("u{0}".format(i), "x[{0}]".format(i + 5)) for i in range(3)] +\
            [("qd{0}".format(i), "dxdt[{0}]".format(i)) for i in range(5)] +\
            [("ud{0}".format(i), "dxdt[{0}]".format(i + 5)) for i in range(3)]

for old, new in state_rep:
  output = output.replace(old, new)

f = file("rattleback_output.txt", 'w')
f.write(output)
f.close()
