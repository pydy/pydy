#!/usr/bin/env python

from ..bodies import Ground

class TestBodies():

    def setup(self):
        self.ground = Ground()

    def test_init(self):
        # Other properties are tested in sympy.RigidBody
