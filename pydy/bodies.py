from sympy.physics.mechanics import RigidBody

from .force import Force

__all__ = ['Body', 'Ground']


class Body(RigidBody):
    """
    Paramters
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
    def __init__(self, name, masscenter=None, frame=None, mass=None, inertia=None):
        self.name = name
        self.parent = None
        self.child = None
        self.force_list = []

        # TODO define properties for RigidBody.

        RigidBody.__init__(self, name, masscenter, frame, mass, inertia)

    def add_force(self, force):
        """
        Adds force to the body by adding Force's instance to the force_list.
        force_list is used by system to get the tuples from the instance and
        create force_list for Kane's method.

        Paramters
        ---------
        force: pydy.Force
            a constant from pydy.Force which will be used by the system to get the exact
            value of the force.

        Example
        -------
        >>> body = Body('body')
        >>> force = Force(body.masscenter, body.frame.z, Symbol('g'))
        >>> body.add_force(force)

        """
        # TODO


class Ground(Body):
    def __init__(self, name, masscenter=None, frame=None, mass=None, inertia=None):
        # TODO change properties of Body.
        Body.__init__(self, name, masscenter, frame, mass, inertia)
