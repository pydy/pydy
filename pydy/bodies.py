from sympy.physics.mechanics import RigidBody, ReferenceFrame, Point, inertia
from sympy import Symbol

__all__ = ['Ground']


class Ground(RigidBody):
    """TODO"""
    def __init__(self):
        name = "Ground"
        reference_frame = ReferenceFrame(name + "_Frame")
        center_of_mass = Point(name + "_CenterOfMass")
        center_of_mass.set_vel(reference_frame, 0)
        ground_inertia = inertia(reference_frame, 1, 1, 1)
        inertia_tuple = (ground_inertia, center_of_mass)
        mass = Symbol(name + "_Mass")

        RigidBody.__init__(self, name, center_of_mass, reference_frame,
                           mass, inertia_tuple)
