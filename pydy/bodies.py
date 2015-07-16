from sympy import Symbol, oo
from sympy.physics.mechanics import RigidBody, Particle, ReferenceFrame, outer
from sympy.physics.vector import Point

__all__ = ['Body', 'Ground']


class Body(object):
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
    inertia : (Dyadic, Point) (optional)
        The body's inertia about a point; stored in a tuple as shown above.

    Example:
    --------
    1. Default behaviour. It creates a RigidBody after defining mass,
     masscenter, frame and inertia.

    >>> from pydy.bodies import Body
    >>> body = Body('name_of_body')

    2. Passing attributes of Rigidbody. All the arguments needed to create a
     RigidBody can be passed while creating a Body too.

    >>> from sympy import Symbol
    >>> from sympy.physics.mechanics import ReferenceFrame, Point
    >>> from pydy.bodies import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Point('masscenter')
    >>> frame = ReferenceFrame('frame')
    >>> I = outer (A.x, A.x)
    >>> inertia_tuple = (I, P)
    >>> body = Body('name_of_body', masscenter, mass, frame, inertia_tuple)

    3. Creating a Particle. If masscenter and mass are passed, and inertia is
     not then a Particle is created.

    >>> from sympy import Symbol
    >>> from sympy import Point
    >>> from pydy.bodies import Body
    >>> mass = Symbol('mass')
    >>> masscenter = Symbol('masscenter')
    >>> body = Body('name_of_body', masscenter, mass)

    Similarly, A frame can also be passed while creating a Particle.

    """
    def __init__(self, name, masscenter=None, mass=None, frame=None,
                 inertia=None):

        self.name = name
        self.parent = None
        self.child = None
        self.force_list = []
        self.body = None  # TODO assign a better name.
        self._coordinates = []
        self._speeds = []
        self._counter = 0

        if masscenter is None:
            self._masscenter = Point(self.name + '_masscenter')
        else:
            self._masscenter = masscenter

        if mass is None:
            self._mass = Symbol(self.name + '_mass')
        else:
            self._mass = mass

        if frame is None:
            self._frame = ReferenceFrame(self.name + '_frame')
        else:
            self._frame = frame

        if inertia is None:
            inertia = outer(self._frame.x, self._frame.x)  # tensor product
            self._inertia = (inertia, self._masscenter)
        else:
            self._inertia = inertia

        # If user passes masscenter and mass then a particle is created
        # otherwise a rigidbody. As a result a body may or may not have inertia.
        if inertia is None and mass is not None:
            self.body = Particle(self.name, self._masscenter, self._mass)
        else:
            self.body = RigidBody(self.name, self._masscenter, self._frame,
                                  self._mass, self._inertia)

    def add_coordinate(self, coordinate):
        self._coordinates.append(coordinate)

    def add_speed(self, speed):
        self._speeds.append(speed)

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
        point_vector = self._convert_tuple_to_vector(point_vector)
        force_vector = self._convert_tuple_to_vector(force_vector)
        point = self._masscenter.locatenew(self.name + '_point' + self._counter, point_vector)
        self.force_list.append((point, force_vector))
        self._counter += 1

    def _convert_tuple_to_vector(self, pos_tuple):
        if len(pos_tuple) == 3:
            unit_vectors = [self.frame.x, self.frame.y, self.frame.z]
            vector = 0
            for i in range(3):
                vector += pos_tuple[i] * unit_vectors[i]
            return vector
        else:
            raise TypeError('position tuple must be of length 3')
