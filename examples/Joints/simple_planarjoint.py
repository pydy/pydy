from sympy import Symbol
from sympy.physics.vector import ReferenceFrame, Point, Vector
from sympy.physics.mechanics import inertia, RigidBody, KanesMethod, \
    dynamicsymbols

# parent
p_name = "parent"
p_frame = ReferenceFrame(p_name + '_frame')
p_masscenter = Point(p_name + '_masscenter')
p_masscenter.set_vel(p_frame, 0)
p_mass = Symbol(p_name + '_mass')
p_inertia = (inertia(p_frame, 1, 1, 1), p_masscenter)
parent = RigidBody(p_name, p_masscenter, p_frame, p_mass, p_inertia)

# child
c_name = "child"
c_frame = ReferenceFrame(c_name + '_frame')
c_masscenter = Point(c_name + '_masscenter')
c_masscenter.set_vel(c_frame, 0)
c_mass = Symbol(c_name + '_mass')
c_inertia = (inertia(c_frame, 1, 1, 1), c_masscenter)
child = RigidBody(c_name, c_masscenter, c_frame, c_mass, c_inertia)

# planar joint
# ------------
# generalized coordinates in specific order.
theta = dynamicsymbols('theta')  # rotation around z axis.
thetad = dynamicsymbols('theta', 1)
omega = dynamicsymbols('omega')
disx = dynamicsymbols('disx')  # translation along x axis.
disxd = dynamicsymbols('disx', 1)
velx = dynamicsymbols('velxx')
disy = dynamicsymbols('disy')  # translation along y axis.
disyd = dynamicsymbols('disy', 1)
vely = dynamicsymbols('vely')

c_frame.orient(p_frame, 'Axis', [0, p_frame.x])

p_joint_point = p_masscenter.locatenew(
    p_name + '_parent_joint',
    Vector(0))

c_joint_point = c_masscenter.locatenew(
    c_name + '_child_joint',
    c_frame.x + c_frame.y + c_frame.z)

c_joint_point.set_pos(p_joint_point, 0)

# Adding rotation
c_frame.orient(p_frame, 'Axis', [theta, p_frame.z])
c_frame.set_ang_vel(p_frame, omega * p_frame.z)

# Adding translation along x axis.
c_joint_point.set_pos(p_joint_point, disx * p_frame.x)
c_joint_point.set_vel(p_frame, velx * p_frame.x)

# Adding translation along y axis
c_joint_point.set_pos(p_joint_point, disy * p_frame.y)
c_joint_point.set_vel(p_frame, vely * p_frame.y)

c_masscenter.v2pt_theory(p_masscenter, p_frame, c_frame)

#Joints Method
# ---------------
q_ind = [theta, disx, disy]
u_ind = [omega, velx, vely]

kd = [thetad - omega, disxd - velx, disyd - vely]
BL = [parent, child]
k = Symbol('k')
FL = [(c_masscenter, k * disx * c_frame.x),
      (c_masscenter, k * disy * c_frame.y)]

KM = KanesMethod(p_frame, q_ind=q_ind, u_ind=u_ind, kd_eqs=kd)
print BL
print FL
print kd
print q_ind
print u_ind
print KM.kanes_equations(FL, BL)
