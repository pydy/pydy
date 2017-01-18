from sympy import Symbol

__all__ = ['Force', 'GravitationalForce']


class Force(object):
    """
    It provides an concept of force, although actual force is added by system.

    Parameters
    ---------
    point: sympy.physics.vector.Point
        Point of application of force.
    direction: sympy.physics.vector.Vector
        Direction of application of force.
    magnitude: (optional)
        defines the magnitude of force.

    Example
    -------
    >>> body = Body('body')
    >>> force = Force(body.masscenter, body.frame.z, Symbol('g'))
    >>> body.add_force(force)

    """
    def __init__(self, point, direction, magnitude=None):
        self.point = point
        self.direction = direction.normalize()
        if magnitude is None:
            self.magnitude = Symbol("ForceMagnitude_" + point.name)
        else:
            self.magnitude = magnitude

    def get_force_tuple(self):
        return (self.magnitude * self.direction, self.point)
