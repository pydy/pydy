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

# spherical joint
thetax = dynamicsymbols('thetax')
thetay = dynamicsymbols('thetay')
thetaz = dynamicsymbols('thetaz')
thetaxd = dynamicsymbols('thetax', 1)
thetayd = dynamicsymbols('thetay', 1)
thetazd = dynamicsymbols('thetaz', 1)
omegax = dynamicsymbols('omegax')
omegay = dynamicsymbols('omegay')
omegaz = dynamicsymbols('omegaz')


p_joint_point = p_masscenter.locatenew(
    p_name + '_parent_joint',
    Vector(0))

c_joint_point = c_masscenter.locatenew(
    c_name + '_child_joint',
    c_frame.x + c_frame.y + c_frame.z)

c_joint_point.set_pos(p_joint_point, 0)

c_frame.orient(p_frame, 'Axis', [thetax, p_frame.x])
c_frame.set_ang_vel(p_frame, omegax * p_frame.x)

c_frame.orient(p_frame, 'Axis', [thetay, p_frame.y])
c_frame.set_ang_vel(p_frame, omegay * p_frame.y)

c_frame.orient(p_frame, 'Axis', [thetaz, p_frame.z])
c_frame.set_ang_vel(p_frame, omegaz * p_frame.z)

c_masscenter.v2pt_theory(p_masscenter, p_frame, c_frame)

# JointsMethod
q_ind = [thetax, thetay, thetaz]
u_ind = [omegax, omegay, omegaz]
kd = [thetaxd - omegax, thetayd - omegay, thetazd - omegaz]
BL = [parent, child]
gravity = Symbol('gravity')
FL = [(c_masscenter, c_mass * gravity * p_frame.y)]

KM = KanesMethod(p_frame, q_ind=q_ind, u_ind=u_ind, kd_eqs=kd)
print BL
print FL
print kd
print q_ind
print u_ind
print KM.kanes_equations(FL, BL)
