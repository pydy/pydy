from sympy import Symbol
from sympy.physics.mechanics import RigidBody, Particle, ReferenceFrame, \
    outer, inertia
from sympy.physics.vector import Point, Vector

__all__ = ['Body']


class Body(RigidBody, Particle):
    """
    A Body which can be connected by joints.

    It can create a Particle or RigidBody depending upon the arguments passed.
    If nothing is passed, then by default, a RigidBody is created. User can
    also pass masscenter, mass, frame and inertia to change the arguments
    passed to RigidBody.
    To create a Particle only mass and masscenter needs to be passed. User can
    also pass a frame which will be attached to the Particle.

    In short, it creates a Rigidbody if inertia is passed or nothing is passed.
    and a Particle if inertia is not passed and things like masscenter and mass
    are.

    Parameters
    ---------
    name: string
        defines the name of the body. It is used as the base for defining body
        specific properties.
    masscenter : Point (optional)
        The point which represents the center of mass of the rigid body.
    frame : ReferenceFrame (optional)
        The ReferenceFrame which the rigid body is fixed in.
    mass : Sympifyable (optional)
        The body's mass.
    body_inertia : Dyadic (instance of inertia)
        The body's inertia about center of mass.

    Example:
    --------
    1. Default behaviour. It creates a RigidBody after defining mass,
     masscenter, frame and inertia.

    >>> from pydy.bodies import Body
    >>> body = Body('name_of_body')

    2. Passing attributes of Rigidbody. All the arguments needed to create a
     RigidBody can be passed while creating a Body too.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import ReferenceFrame, Point, inertia
    >>> from pydy.bodies import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Point('masscenter')
    >>> frame = ReferenceFrame('frame')
    >>> body_inertia = inertia(frame, 1, 0, 0)
    >>> body = Body('name_of_body', masscenter, mass, frame, body_inertia)

    3. Creating a Particle. If masscenter and mass are passed, and inertia is
     not then a Particle is created.

    >>> from sympy import Symbol
    >>> from sympy import Point
    >>> from pydy.bodies import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Point('masscenter')
    >>> body = Body('name_of_body', masscenter, mass)

    Similarly, A frame can also be passed while creating a Particle.

    """
    def __init__(self, name, masscenter=None, mass=None, frame=None,
                 body_inertia=None):

        _name = name
        self.parent = None
        self.child = None
        self.force_list = []
        self._counter = 0

        if masscenter is None:
            self._masscenter = Point(_name + '_masscenter')
        else:
            self._masscenter = masscenter

        if mass is None:
            _mass = Symbol(_name + '_mass')
        else:
            _mass = mass

        if frame is None:
            self._frame = ReferenceFrame(_name + '_frame')
        else:
            self._frame = frame

        if body_inertia is None and mass is None:
            _inertia = (inertia(self._frame, 1, 1, 1), self._masscenter)
        else:
            _inertia = (body_inertia, self._masscenter)

        # If user passes masscenter and mass then a particle is created
        # otherwise a rigidbody. As a result a body may or may not have inertia.
        if body_inertia is None and mass is not None:
            Particle.__init__(self, name, self._masscenter, _mass)
        else:
            RigidBody.__init__(self, _name, self._masscenter, self._frame, _mass, _inertia)

    def add_force(self, point_vector, force_vector):
        """
        Adds force to the body by adding Force's instance to the force_list.
        force_list is used by system to get the tuples from the instance and
        create force_list for Kane's method.

        Parameters
        ---------
        point_vector: A 3 Tuple
            Defines an arbitary point on a body on which the force will be
            applied. The tuple defines the x, y and z values of the point
            w.r.t body's masscenter.
        force_vector: A 3 Tuple
            Defines the force vector w.r.t frame of the body. The tuple defines
            the values of vector in x, y and z directions of body's frame.

        Example
        --------
        To add a force of magnitude 1 in y direction on a point at distance of
        1 in x direction. All the directions are w.r.t body's frame.

        >>> body = Body('body')
        >>> body.add_force((1,0,0), (0,1,0))

        """
        if not isinstance(point_vector, tuple):
            raise TypeError("Point Vector must be a tuple of length 3")
        else:
            point_vector = self._convert_tuple_to_vector(point_vector)

        if not isinstance(force_vector, tuple):
            raise TypeError("Force vector must be a tuple of length 3")
        else:
            force_vector = self._convert_tuple_to_vector(force_vector)

        point = self._masscenter.locatenew(self._name + '_point' + str(self._counter),
                                           point_vector)
        self.force_list.append((point, force_vector))
        self._counter += 1

    def _convert_tuple_to_vector(self, pos_tuple):
        if len(pos_tuple) != 3:
            raise TypeError('position tuple must be of length 3')
        else:
            unit_vectors = [self._frame.x, self._frame.y, self._frame.z]
            vector = Vector(0)
            for i in range(3):
                vector += pos_tuple[i] * unit_vectors[i]
            return vector
