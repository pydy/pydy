from sympy import Symbol
from sympy.physics.vector import ReferenceFrame, Point, Vector
from sympy.physics.mechanics import inertia, RigidBody, KanesMethod, \
    dynamicsymbols

from pydy.bodies import Body

# parent body
#p_name = "parent"
#parent = Body('parent')
#p_frame = parent.get_frame()
#p_masscenter = parent.get_masscenter()
#p_mass = parent.get_mass()

p_name = "parent"
p_frame = ReferenceFrame(p_name + '_frame')
p_masscenter = Point(p_name + '_masscenter')
p_masscenter.set_vel(p_frame, 0)
p_mass = Symbol(p_name + '_mass')
p_inertia = (inertia(p_frame, 1, 1, 1), p_masscenter)
parent = RigidBody(p_name, p_masscenter, p_frame, p_mass, p_inertia)


# child body
# Body can also be used in its place.
#c_name = "child"
#child = Body('child')
#c_frame = child.get_frame()
#c_masscenter = child.get_masscenter()
#c_mass = child.get_mass()

c_name = "child"
c_frame = ReferenceFrame(c_name + '_frame')
c_masscenter = Point(c_name + '_masscenter')
c_masscenter.set_vel(c_frame, 0)
c_mass = Symbol(c_name + '_mass')
c_inertia = (inertia(c_frame, 1, 1, 1), c_masscenter)
child = RigidBody(c_name, c_masscenter, c_frame, c_mass, c_inertia)



# pinjoint
theta = dynamicsymbols('theta')
thetad = dynamicsymbols('theta', 1)
omega = dynamicsymbols('omega')

p_axis = p_frame.x
c_axis = c_frame.x

c_frame.orient(p_frame, 'Axis', [theta, p_axis])
c_frame.set_ang_vel(p_frame, omega * p_axis)

p_joint_point = p_masscenter.locatenew(
    p_name + '_parent_joint',
    Vector(0))

c_joint_point = c_masscenter.locatenew(
    c_name + '_child_joint',
    c_frame.x + c_frame.y + c_frame.z)

c_joint_point.set_pos(p_joint_point, 0)
c_masscenter.v2pt_theory(p_masscenter, p_frame, c_frame)

# JointsMethod
q_ind = [theta]
u_ind = [omega]
kd = [thetad - omega]
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
