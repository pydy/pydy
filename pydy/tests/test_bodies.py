#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
from sympy import Symbol
from sympy.physics.vector import Point, Dyadic, ReferenceFrame
from sympy.physics.mechanics import Particle, RigidBody,\
    dynamicsymbols, outer

# local
from ..bodies import Body


def test_points(point1, point2):
    assert point1.name == point2.name
    assert point1._pos_dict == point2._pos_dict
    assert point1._vel_dict == point2._vel_dict
    assert point1._acc_dict == point2._acc_dict


def test_referenceframes(frame1, frame2):
    assert frame1.name == frame2.name
    assert frame1.str_vecs == frame2.str_vecs
    assert frame1.indices == frame2.indices
    assert frame1.varlist == frame2.varlist
    assert frame1._dlist == frame2._dlist
    assert frame1._ang_vel_dict == frame2._ang_vel_dict
    assert frame1._ang_acc_dict == frame2._ang_acc_dict
    assert frame1._dcm_dict == frame2._dcm_dict

class TestBody():
    def test_default(self):
        body = Body('body')
        assert body.name == 'body'
        assert body.parent is None
        assert body.child is None
        assert body.force_list == []
        assert body._coordinates == []
        assert body._speeds == []
        test_points(body._masscenter, Point('body_masscenter'))
        assert body._mass == Symbol('body_mass')
        test_referenceframes(body._frame, ReferenceFrame('body_frame'))
        assert body._inertia == (outer(body._frame.x, body._frame.x),
                                 body._masscenter)
        assert isinstance(body.body, RigidBody)

        rigidbody = body.body
        test_points(rigidbody.masscenter, body._masscenter)
        assert rigidbody.mass == body._mass
        test_referenceframes(rigidbody.frame, body._frame)
        assert rigidbody.inertia == body._inertia

    def test_custom_rigid_body(self):
        masscenter = Point('rigidbody_masscenter')
        mass = Symbol('rigidbody_mass')
        frame = ReferenceFrame('rigidbody_frame')
        inertia = outer(frame.x, frame.x)
        inertia_tuple = (inertia, masscenter)
        self.rigidbody_body = Body('rigidbody_body', masscenter, mass, frame,
                                   inertia_tuple)
        test_points(self.rigidbody_body._masscenter, masscenter)
        assert self.rigidbody_body._mass == mass
        test_referenceframes(self.rigidbody_body._frame, frame)
        assert self.rigidbody_body._inertia == inertia_tuple
        assert isinstance(self.rigidbody_body.body, RigidBody)

    def test_particle_body(self):
        masscenter = Point('particle_masscenter')
        mass = Symbol('particle_mass')
        frame = ReferenceFrame('particle_frame')
        self.particle_body = Body('particle_body', masscenter, mass, frame)
        test_points(self.particle_body._masscenter, masscenter)
        assert self.particle_body._mass == mass
        test_referenceframes(self.particle_body._frame, frame)
        assert self.particle_body._inertia is None
        assert isinstance(self.particle_body.body, Particle)

    def test_add_coordinates(self):
        for body in [self.rigidbody_body, self.particle_body]:
            assert body._coordinates == []
            q1 = dynamicsymbols('q1')
            body.add_coordinate(q1)
            assert q1 in body._coordinates

    def test_add_speeds(self):
        for body in [self.rigidbody_body, self.particle_body]:
            assert body._speeds == []
            u1 = dynamicsymbols('u1')
            body.add_speed(u1)
            assert u1 in body._speeds

    def test_particle_body_add_force(self):
        a = Symbol('a')
        self.particle_body.add_force((0,0,0), (a,0,0))
        assert len(self.particle_body.force_list) == 1
        point = self.particle_body._masscenter.locatenew(
            self.particle_body.name + '_point0', 0)
        force_vector = a * self.particle_body._frame.x
        test_points(self.particle_body.force_list[0][0], point)
        assert self.particle_body.force_list[0][1] == force_vector

    def test_body_add_force(self):
        l = Symbol('l')
        Fa = Symbol('Fa')
        point = self.rigidbody_body._masscenter.locatenew(
            'rigidbody_body_point0',
            l * self.rigidbody_body._frame.x)
        self.rigidbody_body.add_force((l,0,0), (0,0,Fa))
        assert len(self.rigidbody_body.force_list) == 1
        force_vector = Fa * self.rigidbody_body._frame.z
        test_points(self.rigidbody_body.force_list[0][0], point)
        assert self.rigidbody_body.force_list[0][1] == force_vector
