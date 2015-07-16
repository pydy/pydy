#!/usr/bin/env python
# -*- coding: utf-8 -*-

# external libraries
from sympy import Symbol
from sympy.physics.vector import Point, Dyadic
from sympy.physics.mechanics import ReferenceFrame, Particle, RigidBody, dynamicsymbols

# local
from ..bodies import Body, Ground


class TestBody():
    def setup(self):
        self.body = Body('Body')

    def test_body_init(self):
        # by default a rigidbody is created.
        assert self.body.name == 'Body'
        assert self.body.parent is None
        assert self.body.child is None
        assert isinstance(self.body._masscenter, Point)
        assert isinstance(self.body._frame,    ReferenceFrame)
        assert isinstance(self.body._mass, Symbol)
        assert len(self.body._inertia) == 2
        assert isinstance(self.body._inertia[0], Dyadic)
        assert isinstance(self.body._inertia[1], Point)
        assert isinstance(self.body.body, RigidBody)

    def test_body_properties_value(self):
        name = self.body.name
        assert self.body._masscenter == Point(name + ' MassCenter')
        assert self.body._frame == ReferenceFrame(name + ' ReferenceFrame')
        assert self.body._mass == Symbol(name + ' Mass')

    def test_body_particle(self):
        point = Point('point')
        mass = Symbol('mass')
        self.particle_body = Body(masscenter=point, mass=mass)
        assert isinstance(self.particle_body.body, Particle)
        assert self.particle_body._masscenter == point
        assert self.particle_body._mass == mass

    def test_body_particle_add_force(self):
        vector = Symbol('a') * self.particle_body._frame.x
        self.particle_body.add_force(vector)
        assert self.particle_body.force_list == [(self.particle_body._masscenter, vector)]

    def test_body_add_force(self):
        assert isinstance(self.body.force_list, list())
        assert len(self.body.force_list) == 0

        point = self.body._masscenter.locatenew('force_point', Symbol('l')*self.body._frame.x)
        force_tuple = (point, Symbol('Fa')*self.body._frame.z)
        self.body.add_force(force_tuple)
        assert len(self.body.force_list) == 1
        assert self.body.force_list[0] == [force_tuple]

    def test_body_coordinates(self):
        assert self.body._coordinates == []
        q1 = dynamicsymbols('q1')
        self.body.add_coordinate(q1)
        assert q1 in self.body._coordinates

    def test_body_speeds(self):
        assert self.body._speeds == []
        u1 = dynamicsymbols('u1')
        self.body.add_speed(u1)
        assert u1 in self.body._speeds


class TestGround():
    def setup(self):
        self.ground = Ground()

    def test_ground_init(self):
        assert hasattr(self.ground, 'origin')
        assert self.ground._masscenter == Point(self.ground.name + '_Origin')
        assert self.ground.parent is None
        assert self.ground._frame == ReferenceFrame(self.ground.name + '_Frame')
