from sympy import Symbol
from sympy.physics.vector import Point
from sympy.physics.mechanics import dynamicsymbols

from ..bodies import Body
from ..joints import Joint


class TestJoints():
    def setup(self):
        self.parent = Body('parent')
        self.child = Body('child')
        self.joint = Joint('joint', self.parent, self.parent.frame.x, self.child, 0)

    def test_joint_init(self):
        assert self.joint.name == 'joint'
        assert self.joint.child == self.child
        assert self.joint.parent == self.parent
        assert self.joint.parent_joint_vector == self.parent.frame.x
        assert self.joint.child_joint_vector == 0

    def test_joint_set_parent_child_rel(self):
        self.joint.set_parent_child_rel()
        assert self.child.parent == self.parent
        assert self.parent.child == self.child

    def test_joint_locate_joint_point(self):
        self.joint._locate_joint_point()
        assert hasattr(self.joint, 'joint_point')

        joint_point = self.parent.masscenter.locatenew(self.name + '_Point',
                                                       self.joint.parent_joint_vector)
        assert self.joint_point.name == joint_point

        child_masscenter = Point(self.child.name + '_MassCenter')
        child_masscenter.set_pos(joint_point, self.joint.child_joint_vector)
        assert self.child.masscenter == child_masscenter

    def test_joint_set_child_vel(self):
        # Note: this numeric and symbolic value implementation can be applied anywhere in joints
        # wherever symbol is used e.g. in set_joint_point_vel()

        # 1. numeric value
        self.joint.set_child_vel(10 * self.parent.frame.x)
        assert self.child.masscenter.vel(self.parent.frame) == 10 * self.parent.frame.x

        # 2. symbolic value
        v1 = dynamicsymbols('v1')
        self.joint.set_child_vel(v1 * self.parent.frame.x)
        assert self.child.masscenter.vel(self.parent.frame) == v1 * self.parent.frame.x

    def test_joint_set_joint_point_vel(self):
        v2 = dynamicsymbols('v2')
        self.joint.set_joint_point_vel(v2 * self.parent.frame.x)
        assert self.joint.joint_point.vel(self.parent.frame) == v2 * self.parennt.frame.x

    def test_joint_set_child_angular_velocity(self):
        v3 = dynamicsymbols('v3')
        self.joint.set_child_angular_vel(v3 * self.parent.frame.z)
        assert self.child.frame.ang_vel_in(self.parent.frame) == v3 * self.parent.frame.z

    def test_joint_orient_child(self):
        q1 = Symbol('q1')
        self.joint.orient_child('Axis', [q1, self.parent.frame.x])

        # TODO other similar frame.orient() tests can be done.
        # leaving more tests to discuss its necessity with others.


class TestPinJoint():
    def setup(self):
        #TODO
