from sympy.physics.mechanics import RigidBody

__all__ = ['Body', 'Ground']

class Body(RigidBody):
    """Base class of every body"""
    def __init__(self, name):
        self.name = name
        self.reference_frame = ReferenceFrame(self.name + "_Frame")
        self.center_of_mass = Point(self.name + "_CenterOfMass")
        self.center_of_mass.set_vel(self.reference_frame, 0)
        self.ground_inertia = inertia(self.reference_frame, 1, 1, 1)
        self.inertia_tuple = (self.ground_inertia, self.center_of_mass)
        self.mass = Symbol(self.name + "_Mass")

        RigidBody.__init__(self, self.name, self.center_of_mass,
                           self.reference_frame, self.mass,
                           self.inertia_tuple)

        self.parent = None
        self.child = None

class Ground(Body):
    """Root body in every pydy.system.System when created using joints.

    """
    def __init__(self, name):
        Body.__init__(self, name)
        # TODO
