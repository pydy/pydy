from sympy import Symbol
from sympy.physics.vector import ReferenceFrame, Point
from sympy.physics.mechanics import inertia, RigidBody, KanesMethod, \
    dynamicsymbols

from pydy.bodies import Body
from pydy.joints import PinJoint
from pydy.joints_method import JointsMethod


class TestPinJoint():
    def setup(self):
        parent = Body('parent')
        child = Body('child')

        gravity = Symbol('gravity')
        # add gravity force in y direction at the center of mass of child.
        child.add_force((0, child.get_mass() * gravity, 0), child.get_masscenter())

        self.pin_joint = PinJoint('pin_joint', parent, child, parent_point_pos=(1, 1, 1),
                             child_point_pos=(1, 1, 1))

        JM = JointsMethod([self.pin_joint], parent)
        (self.fr, self.frstar) = JM.get_kanes()

    def test_equations_of_motion(self):
        # parent
        p_frame = ReferenceFrame(p_name + '_frame')
        p_masscenter = Point(p_name + '_masscenter')
        p_masscenter.set_vel(p_frame, 0)
        p_mass = Symbol(p_name + '_mass')
        p_inertia = (inertia(p_frame, 1, 1, 1), p_masscenter)
        parent = RigidBody(p_name, p_masscenter, p_frame, p_mass, p_inertia)
        # child
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
            p_frame.x + p_frame.y + p_frame.z)

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
        FL = [(c_masscenter, c_mass * gravity * c_frame.y)]

        KM = KanesMethod(p_frame, q_ind=q_ind, u_ind=u_ind, kd_eqs=kd)
        (fr, frstar) = KM.kanes_equations(FL, BL)

        assert self.fr == fr
        assert self.frstar == frstar
