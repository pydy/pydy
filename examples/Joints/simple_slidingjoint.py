from sympy import Symbol
from sympy.physics.vector import ReferenceFrame, Point
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

# sliding joint
dis = dynamicsymbols('dis')
disd = dynamicsymbols('dis', 1)
vel = dynamicsymbols('vel')

p_direction = p_frame.x
c_direction = c_frame.x

c_frame.orient(p_frame, 'Axis', [0, p_frame.x])

p_joint_point = p_masscenter.locatenew(
    p_name + '_parent_joint',
    p_frame.x + p_frame.y + p_frame.z)

c_joint_point = c_masscenter.locatenew(
    c_name + '_child_joint',
    c_frame.x + c_frame.y + c_frame.z)

p_joint_point.set_vel(p_frame, 0)
c_joint_point.set_vel(c_frame, 0)

c_joint_point.set_pos(p_joint_point, dis * p_direction)
c_joint_point.set_vel(p_frame, vel * p_direction)
c_masscenter.set_vel(p_frame, vel * p_direction)

# JointsMethod
q_ind = [dis]
u_ind = [vel]
kd = [disd - vel]
BL = [parent, child]
gravity = Symbol('gravity')
spring_constant = Symbol('k')
# Note: Here, spring force is in direction c_frame.x + c_frame.y
# but due to sliding joint's restriction of degree of freedom,
# the child is moving only in c_frame.x
FL = [(c_masscenter, c_mass * gravity * c_frame.x),
      (c_masscenter, spring_constant * dis * (c_frame.x + c_frame.y))]
KM = KanesMethod(p_frame, q_ind=q_ind, u_ind=u_ind, kd_eqs=kd)
print BL
print FL
print kd
print q_ind
print u_ind
print KM.kanes_equations(FL, BL)
