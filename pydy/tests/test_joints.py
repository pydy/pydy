from sympy import Symbol, acos
from sympy.physics.vector import Point, dot, cross
from sympy.physics.mechanics import dynamicsymbols

from ..bodies import Body
from ..joints import Joint, PinJoint, SlidingJoint


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
        self.joint._set_parent_child_rel()
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
        self.parent = Body('parent')
        self.child = Body('child')
        self.pinjoint = PinJoint('pinjoint', self.parent, self.child)

    def test_pinjoint_init(self):
        # default values
        assert self.pinjoint.parent_axis == self.parent.frame.x
        assert self.pinjoint.child_axis == self.child.frame.x
        assert self.pinjoint.parent_point == self.parent.masscenter
        assert self.pinjoint.child_point == self.child.masscenter
        assert self.pinjoint.parent_joint_vector == self.parent.frame.x
        assert self.pinjoint.child_joint_vector == self.child.frame.x

        assert self.parent_point.vel(self.parent.frame) == 0
        assert self.child_point.vel(self.child.frame) == 0

    def test_pinjoint_init_using_particles(self):
        point1 = Point('point1')
        point2 = Point('point2')
        mass1 = Symbol('m1')
        mass2 = Symbol('m2')

        parent = Particle('parent', point1, mass1)
        child = Particle('child', point2, mass2)

        pinjoint = PinJoint('pinjoint', parent, child)
        assert pinjoint.parent = parent
        assert pinjoint.child = child
        assert

    def test_pinjoint_parameters(self):
        self.pinjoint = PinJoint('pinjoint', self.parent, self.child,
                                 (1, 0, 0), (0, 1, 0), parent_axis='x', child_axis='y')
        assert self.pinjoint.parent_axis == self.parent.frame.x
        assert self.pinjoint.child_axis == self.child.frame.y
        point1 = self.parent.masscenter.locatenew(self.parent.frame, self.parent.frame.x)
        point2 = self.child.masscenter.locatenew(self.child.frame, self.child.frame.y)
        assert self.pinjoint.parent_point == point1
        assert self.pinjoint.child_joint_piont == point2
        assert self.pinjoint.parent_joint_vector == self.parent.frame.x
        assert self.pinjoint.child_joint_vector == self.chidl.frame.y

    def test_pinjoint_tuple_vector(self):
        a, b, c, d, e, f = symbols('a b c d e f')
        pinjoint = PinJoint('pinjoint', self.parent, self.child,
                            (a, b, c), (d, e, f), parent_axis='x', child_axis='y')



    def test_pinjoint_functions(self):

        # part 1 assining parent-child relationship
        self.pinjoint._set_parent_child_rel()
        assert self.child.parent == self.parent
        assert self.parent.child == self.child

        # part 2 locating joint point in bodies
        self.pinjoint._locate_joint_point()
        child_masscenter = self.child.masscenter
        parent_masscenter = self.parent.masscenter
        parent_point = parent_masscenter.locatenew(self.name + '_Point', self.parent_joint_vector)
        child_point = child_masscenter.locatenew(self.name + '_Point', self.child_joint_vector)
        assert self.pinjoint.parent_joint_vector == parent_point
        assert self.pinjoint.child_joint_vector == child_point

        # part 3 aligning parent_axis and child_axis
        self.join_frames()
        # sample implementation
        # self.child.frame.orient(self.parent.frame, 'Axis', [0, self.parent.frame.x])
        assert cross(self.pinjoint.parent_axis, self.pinjoint.child_axis) == self.child.frame.z
        assert cross(self.pinjoint.child_axis, self.pinjoint.parent_axis) == - self.parent.frame.z
        self.align_axes()
        parent_axis = self.pinjoint.parent_axis
        child_axis = self.pinjoint.child_axis
        # sample implementation
        # perpendicular_axis_in_parent = - cross(parent_axis, child_axis)
        # angle_between_axes = acos(dot(self.pinjoint.child_axis, self.pinjoint.parent_axis))
        # self.child.frame.orient(self.parent.frame, 'Axis',
        #                        [acos(dot(parent_axis, child_axis)/(parent_axis.magnitude() * child_axis.magnitude())),
        #                         cross(parent_axis, child_axis)])
        assert acos(dot(parent_axis, child_axis)/(parent_axis.magnitude() * child_axis.magnitude())) == 0

        # part 4 adding angular velocity to child w.r.t parent
        assert self.child.frame.ang_vel_in(self.parent.frame) != 0

    def test_pinjoint_apply_joint(self):
        # apply_joint() should do everything done above by calling specific funtions
        self.new_pinjoint = PinJoint('pinjoint', self.parent, self.child)
        self.new_pinjoint._apply_joint()
        # part 1
        assert self.child.parent == self.parent
        assert self.parent.child == self.child

        # part 2
        child_masscenter = self.child.masscenter
        parent_masscenter = self.parent.masscenter
        parent_point = parent_masscenter.locatenew(self.name + '_Point', self.parent_joint_vector)
        child_point = child_masscenter.locatenew(self.name + '_Point', self.child_joint_vector)
        assert self.pinjoint.parent_joint_vector == parent_point
        assert self.pinjoint.child_joint_vector == child_point

        # part 4
        parent_axis = self.pinjoint.parent_axis
        child_axis = self.pinjoint.child_axis
        assert acos(dot(parent_axis, child_axis)/(parent_axis.magnitude() * child_axis.magnitude())) == 0


class TestSlidingJoint():
    def setup(self):
        self.parent = Body('parent')
        self.child = Body('child')
        self.slidingjoint = SlidingJoint('slidingjoint', self.parent, self.child)

    def test_slidingjoint_init(self):
        # Default values
        assert self.slidingjoint.parent == self.parent
        assert self.slidingjoint.child == self.child
        assert self.slidingjoint.direction1 == self.parent.x
        assert self.slidingjoint.direction2 == self.child.x
        assert self.slidingjoint.parent_point == self.parent.masscenter
        assert self.slidingjoint.child_point == self.child.masscenter

    def test_slidingjoint_paramters(self):
        # custom parameters
        self.slidingjoint = SlidingJoint('slidingjoint', self.parent, self.child,
                                         (1,0,0), (0,1,0), direction1='x', direction2='y')
        assert self.pinjoint.direction1 == self.parent.frame.x
        assert self.pinjoint.direction2 == self.child.frame.y
        point1 = self.parent.masscenter.locatenew(self.parent.frame, self.parent.frame.x)
        point2 = self.child.masscenter.locatenew(self.child.frame, self.child.frame.y)
        assert self.pinjoint.parent_point == point1
        assert self.pinjoint.child_joint_piont == point2
        assert self.pinjoint.parent_joint_vector == self.parent.frame.x
        assert self.pinjoint.child_joint_vector == self.child.frame.y

    def test_slidingjoint_functions(self):
        self.slidingjoint._set_parent_child_rel()
        assert self.parent.child == self.child
        assert self.child.parent == self.parent

        self.slidingjoint._locate_joint_point()
        child_masscenter = self.child.masscenter
        parent_masscenter = self.parent.masscenter
        parent_point = parent_masscenter.locatenew(self.name + '_Point', self.parent_joint_vector)
        child_point = child_masscenter.locatenew(self.name + '_Point', self.child_joint_vector)
        assert self.pinjoint.parent_joint_vector == parent_point
        assert self.pinjoint.child_joint_vector == child_point

        self.join_frames()
        assert cross(self.pinjoint.parent_axis, self.pinjoint.child_axis) == self.child.frame.z
        assert cross(self.pinjoint.child_axis, self.pinjoint.parent_axis) == - self.parent.frame.z

        self.align_axes()
        parent_axis = self.pinjoint.parent_axis
        child_axis = self.pinjoint.child_axis
        assert acos(dot(parent_axis, child_axis)/(parent_axis.magnitude() * child_axis.magnitude())) == 0

        assert self.child_point.vel(self.parent.frame) != 0


class TestCylindricalJoint():
    def setup(self):
