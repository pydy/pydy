#!/usr/bin/env python

SYMPY_VERSION = sm.__version__

from ..bodies import Body
from ..joints import PinJoint
from ..system import System


class TestPinJoint():
    def setup(self):
        self.sys = System(method='joints')

    def test_joints_parent_child_relationship(self):
        body = Body('body')
        pin_joint = PinJoint(self.ground, 0, body, sys.reference_frame.x * Symbol('l'))

        assert self.sys.ground.child is None
        assert body.parent is None

        self.sys.add_joint(pin_joint)

        assert self.sys.ground.child == body
        assert body.parent == self.sys.ground

    def test_joints_system_body_list(self):
        body = Body('body')
        pin_joint = PinJoint(self.ground, 0, body, sys.reference_frame.x * Symbol('l'))

        assert self.sys.body_list == list()
        self.sys.add_joint(pin_joint)
        assert body in self.sys.body_list

    def test_joint_locate_pin_point(self):
        # TODO

