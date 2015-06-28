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
        Appends the constant from Force class to the force_list class. Applying exact
        value of force is handled by the system.

        Paramters
        ---------
        force: pydy.Force
            a constant from pydy.Force which will be used by the system to get the exact
            value of the force.

        Example
        -------
        >>> body = Body('body')
        >>> body.add_force(Force.GRAVITY)

        """
        # TODO


class Ground(Body):
    def __init__(self, name, masscenter=None, frame=None, mass=None, inertia=None):
        # TODO change properties of Body.
        Body.__init__(self, name, masscenter, frame, mass, inertia)
