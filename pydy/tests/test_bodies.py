#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
from sympy import Symbol
from sympy.physics.vector import Point, Dyadic
from sympy.physics.mechanics import ReferenceFrame, Particle, RigidBody,\
    dynamicsymbols, outer

# local
from ..bodies import Body


class TestBody():
    def test_default(self):
        body = Body('body')
        assert body.name = 'body'
        assert body.parent is None
        assert body.child is None
        assert body.force_list == []
        assert body._coordinates == []
        assert body._speeds == []
        assert body._masscenter == Point('body_masscenter')
        assert body._mass == Symbol('body_mass')
        assert body._frame == ReferenceFrame('body_frame')
        assert body._inertia == (outer(body._frame.x, body._frame.x),
                                 body._masscenter)
        assert isinstance(body.body, RigidBody)

        rigidbody = body.body
        assert rigidbody.masscenter == body._masscenter
        assert rigidbody.mass == body._mass
        assert rigidbody.frame == body._frame
        assert rigidbody.inertia == body._inertia

    def test_custom_rigid_body(self):
        masscenter = Point('rigidbody_masscenter')
        mass = Symbol('rigidbody_mass')
        frame = ReferenceFrame('rigidbody_frame')
        inertia = outer(frame.x, frame.y)
        self.rigidbody_body = Body('rigidbody_body', masscenter, mass, frame, inertia)
        assert self.rigidbody_body._masscenter == masscenter
        assert self.rigidbody_body._mass == mass
        assert self.rigidbody_body._frame == frame
        assert self.rigidbody_body._inertia == inertia
        assert isinstance(self.rigidbody_body.body, RigidBody)

    def test_particle_body(self):
        masscenter = Point('particle_masscenter')
        mass = Symbol('particle_mass')
        frame = ReferenceFrame('particle_frame')
        self.particle_body = Body('particle_body', masscenter, mass, frame)
        assert self.particle_body._masscenter == masscenter
        assert self.particle_body._mass == mass
        assert self.particle_body._frame == frame
        assert self.particle_body._inertia is None
        assert isinstance(self.particle_body.body, Particle)

    def test_add_coordinates(self):
        assert self.body._coordinates == []
        q1 = dynamicsymbols('q1')
        self.body.add_coordinate(q1)
        assert q1 in self.body._coordinates

    def test_add_speeds(self):
        assert self.body._speeds == []
        u1 = dynamicsymbols('u1')
        self.body.add_speed(u1)
        assert u1 in self.body._speeds

    def test_particle_body_add_force(self):
        a = Symbol('a')
        self.particle_body.add_force((0,0,0), (a,0,0))
        assert len(self.particle_body.force_list) == 1
        point = self.particle_body._masscenter
        force_vector = a * self.particle_body._frame.x
        assert self.particle_body.force_list == [(point, force_vector)]

    def test_body_add_force(self):
        l = Symbol('l')
        Fa = Symbol('Fa')
        point = self.rigidbody_body._masscenter.locatenew(
            'rigidbody_body_point0',
            l * self.rigidbody_body._frame.x)
        self.rigidbody_body.add_force((l,0,0), (0,0,Fa))
        assert len(self.rigidbody_body.force_list) == 1
        force_vector = Fa * self.rigidbody_body._frame.z
        assert self.rigidbody_body.force_list == [(point, force_vector)]
