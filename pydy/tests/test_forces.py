from sympy import Symbol
from sympy.physics.vector import Point
from sympy.physics.mechanics import ReferenceFrame

from ..forces import Force

class TestForce():
    def setup(self):
        self.point = Point('point')
        self.reference_frame = ReferenceFrame('referenceframe')
        self.point.set_vel(self.reference_frame, 0)
        self.force = Force(self.point, self.reference_frame.z, Symbol('g'))

    def test_force_init(self):
        assert self.force.point == self.point
        assert self.force.direction == self.direction.normalize()
        assert self.force.magnitude == self.magnitude

    def test_force_get_force_tuple(self):
        force_tuple = self.force.get_force_tuple()
        assert force_tuple = (self.direction.normalize() * self.magnitude, self.point)
