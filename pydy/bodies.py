from sympy.physics.mechanics import RigidBody, Particle

from .force import Force

__all__ = ['Body', 'Ground']


class Body(object):
    """

    Parameters
    ---------
    name: string
        defines the name of the body. It is used as the base for defining body specific
        properties.
    (Rest are same as RigidBody but optional)
    masscenter : Point (optional)
        The point which represents the center of mass of the rigid body.
    frame : ReferenceFrame (optional)
        The ReferenceFrame which the rigid body is fixed in.
    mass : Sympifyable (optional)
        The body's mass.
    inertia : (Dyadic, Point) (optional)
        The body's inertia about a point; stored in a tuple as shown above.

    """
    def __init__(self, name, masscenter=None, mass=None, frame=None, inertia=None):

        self.name = name
        self.parent = None
        self.child = None
        self.force_list = []
        self.body = None # TODO assign a better name.
        self.coordinates = []
        self.speeds = []

        # TODO define properties for RigidBody.

        # is user passes masscenter and mass then a particle is created otherwise a rigidbody.
        # as a result a body may or may not have inertia.
        if inertia is None and mass is not None:
                self.body = Particle(self.name, self.masscenter, self.mass)
        else:
                self.body = RigidBody(self.name, self.masscenter, self.frame, self.mas

    def add_force(self, force):
        """
        Adds force to the body by adding Force's instance to the force_list.
        force_list is used by system to get the tuples from the instance and
        create force_list for Kane's method.

        Paramters
        ---------
        force: (point, vector) or vector
            Adds the vector force to the point. point must have to be in the body.
            If body is RigidBody, pass (point, vector) and if body is Particle then pass vector,
            force is always applied to the masscenter for particle.
        """
        # TODO

    def add_coordinate(self, coordinate):
        # TODO
        self.coordinates.append(coordinate)

    def add_speed(self, speed):
        # TODO
        self.speeds.append(speed)


class Ground(Body):
    def __init__(self):
        # TODO create defualt attributes for Ground.
        super(Ground, self).__init__(*args, **kwargs)
