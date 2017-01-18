#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
from sympy import Symbol
from sympy.physics.vector import Point, ReferenceFrame
from sympy.physics.mechanics import Particle, RigidBody, inertia


# local
from ..bodies import Body

def check_points(point1, point2):
    assert point1.name == point2.name
    assert point1._pos_dict == point2._pos_dict
    assert point1._vel_dict == point2._vel_dict
    assert point1._acc_dict == point2._acc_dict

def check_reference_frames(frame1, frame2):
    assert frame1.name == frame2.name
    assert frame1.str_vecs == frame2.str_vecs
    assert frame1.indices == frame2.indices
    assert frame1.varlist == frame2.varlist
    assert frame1._dlist == frame2._dlist
    assert frame1._ang_vel_dict == frame2._ang_vel_dict
    assert frame1._ang_acc_dict == frame2._ang_acc_dict
    assert frame1._dcm_dict == frame2._dcm_dict


class TestBody():

    def setup(self):
        # Body with RigidBody.
        self.rigidbody_masscenter = Point('rigidbody_masscenter')
        self.rigidbody_mass = Symbol('rigidbody_mass')
        self.rigidbody_frame = ReferenceFrame('rigidbody_frame')
        self.body_inertia = inertia(self.rigidbody_frame, 1, 0, 0)
        self.rigid_body = Body('rigidbody_body', self.rigidbody_masscenter, self.rigidbody_mass,
                               self.rigidbody_frame, self.body_inertia)
        #  Body with Particle
        self.particle_masscenter = Point('particle_masscenter')
        self.particle_mass = Symbol('particle_mass')
        self.particle_frame = ReferenceFrame('particle_frame')
        self.particle_body = Body('particle_body', self.particle_masscenter, self.particle_mass,
                                  self.particle_frame)

    def test_default(self):
        body = Body('body')
        assert body._name == 'body'
        assert body.parent is None
        assert body.child is None
        assert body.force_list == []
        check_points(body.get_masscenter(), Point('body_masscenter'))
        assert body.get_mass() == Symbol('body_mass')
        check_reference_frames(body.get_frame(), ReferenceFrame('body_frame'))
        assert body.get_inertia() == (inertia(body.get_frame(), 1, 1, 1), body.get_masscenter())

    def test_custom_rigid_body(self):
        check_points(self.rigid_body.get_masscenter(), self.rigidbody_masscenter)
        assert self.rigid_body.get_mass() == self.rigidbody_mass
        check_reference_frames(self.rigid_body.get_frame(), self.rigidbody_frame)
        assert self.rigid_body.get_inertia() == (self.body_inertia, self.rigidbody_masscenter)

    def test_particle_body(self):
        check_points(self.particle_body.get_masscenter(), self.particle_masscenter)
        assert self.particle_body.get_mass() == self.particle_mass
        check_reference_frames(self.particle_body.get_frame(), self.particle_frame)
        assert not hasattr(self.particle_body, "_inertia")

    def test_particle_body_add_force(self):
        a = Symbol('a')
        self.particle_body.add_force((0, 0, 0), (a, 0, 0))
        assert len(self.particle_body.force_list) == 1
        point = self.particle_body.get_masscenter().locatenew(
            self.particle_body._name + '_point0', 0)
        force_vector = a * self.particle_body.get_frame().x
        check_points(self.particle_body.force_list[0][0], point)
        assert self.particle_body.force_list[0][1] == force_vector

    def test_body_add_force(self):
        l = Symbol('l')
        Fa = Symbol('Fa')
        point = self.rigid_body.get_masscenter().locatenew(
            'rigidbody_body_point0',
            l * self.rigid_body.get_frame().x)
        self.rigid_body.add_force((l, 0, 0), (0, 0, Fa))
        assert len(self.rigid_body.force_list) == 1
        force_vector = Fa * self.rigid_body.get_frame().z
        check_points(self.rigid_body.force_list[0][0], point)
        assert self.rigid_body.force_list[0][1] == force_vector
