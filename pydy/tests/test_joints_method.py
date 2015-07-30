from ..joints_method import JointsMethod
from ..bodies import Ground, Body
from ..joints import PinJoint, SlidingJoint, CylindricalJoint, SphericalJoint, PlanarJoint


class TestJointsMethod():
    def setup(self):
        self.ground = Ground()
        self.body_1 = Body('first_body')
        self.pinjoint_1 = PinJoint('pinjoint', self.ground, self.body_1, (1,0,0), (0,1,0), parent_axis='x', child_axis='y')
        self.joints_method = JointsMethod(self.ground)

    def test_adding_joint_after_instantiating_jointsmethod(self):
        # you can add joints even after instantiating JointsMethod.
        self.body_2 = Body('second_body')
        self.slidingjoint_1 = SlidingJoint('slidingjoint', self.body_1, self.body_2, (1,-1,0))
        self.joints_method.get_all_bodies()
        assert self.joints_method.bodies == [self.ground, self.body_1, self.body_2]

    def test_get_equations(self):
        self.body_3 = Body('third_body')
        self.pinjoint_2 = PinJoint('pinjoint_1', self.body_2, self.body_3, (1,0,0))
        self.equation = self.joints_method.get_equations()

        assert self.joints_method.bodies == [self.ground, self.body_1, self.body_2, self.body_3]
        # TODO generate actual equations and test
