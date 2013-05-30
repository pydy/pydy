# Sliding Ladder in a Frictionless World

# A ladder is placed on a floor and is leaning against a wall.
# Both of the surfaces in contact with the ladder are frictionless.
# There are two reaction forces on the ladder (one from the floor and
# one from the wall). The ladder falls under the influence of gravity.

# We need only one angle to fully define the position of the ladder.
# The mass, length and moment of intertia of the ladder are given.

# First we define all the symbols we need
# The angle and angular velocity of the ladder
q, u = dynamicsymbols('q u')

# The time derivatives of the symbols declared above
qd, ud = dynamicsymbols('q u', 1)

# The mass, length, and moment of inertia of the ladder and
# the acceleration due to gravity
m, l, g, Izz = symbols('m l g Izz')

# We define the inertial frame and a lean frame for the ladder
N = ReferenceFrame('N')
L = N.orientnew('L', 'Axis', [q, N.z])

# and then set the angular velocity of the lean frame in the inertial
# frame. The angular acceleration is automatically computed
L.set_ang_vel(N, u * N.z)

# Now the origin and the center of mass of the ladder are defined
O = Point('O')
A = Point('A')

# and we use the length and angle to locate the center of mass relative
# to the origin
A.set_pos(O, -l / 2 * cos(q) * N.x + l / 2 * sin(q) * N.y)

# Take the derivative of the position of the center of mass to get the
# corresponding velocity
O.set_vel(N, 0)
A.set_vel(N, l / 2 * u * sin(q) * N.x + l / 2 * u * cos(q) * N.y)

# The ladder can now be defined as a rigid body
ladder = RigidBody('ladder', A, L, m, (inertia(L, 0, 0, Izz), A))

# Set up all the inputs to Kanes Method
kd = [u - qd]
bodyList = [ladder]
forceList = [(A, -m * g * N.y)]

# Finally we solve the dynamics
KM = KanesMethod(N, q_ind = [q], u_ind = [u], kd_eqs = kd)
KM.kanes_equations(forceList, bodyList)

# The mass matrix and the forcing function can be taken out of the
# Kanes Method object KM
MM = KM.mass_matrix
forcing = KM.forcing

# and those can be used to find the equations of motion
rhs = MM.inv() * forcing
kdd = KM.kindiffdict()

rhs = rhs.subs(kdd)
rhs.simplify()
