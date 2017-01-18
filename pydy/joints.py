from sympy import acos, sympify
from sympy.physics.vector import cross, dot, Vector
from sympy.physics.mechanics import dynamicsymbols

__all__ = ['Joint', 'PinJoint', 'SlidingJoint', 'CylindricalJoint',
           'SphericalJoint', 'PlanarJoint']


class Joint(object):
    """Base class for all Joints

    Joints take two bodies as arguments (one is parent body and other is child
    body) and define movements of child w.r.t parent. The movement defined is
    specific to type of Joint. This base class holds all the common and public
    methods to interact with all Joints. User can also create a custom joint
    by creating its subclass and overriding 'apply_joint' method. In this method
    define the movements of child w.r.t to parent.

    Parameters
    ----------
    name: String
        It defines the name of a joint which makes it unique.
    parent: Body
        Parent body.
    child: Body
        Child body.
    parent_point_pos: 3 Tuple (optional)
        Defines the vector to the Joint's point where the parent will be
        connected to child. 3 Tuple defines the values of x, y and z directions
        w.r.t parent's frame.
    child_point_pos: 3 Tuple (optional)
        Defines the vector to the Joint's point where the child will be
        connected to parent. 3 Tuple defines the values of x,y and z directions
        w.r.t child's frame.

    """

    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None):
        self.name = name
        self.parent = parent
        self.child = child
        self._coordinates = []
        self._speed = []
        self._kd = []

        if parent_point_pos is None:
            parent_point_pos = (0, 0, 0)  # Center of mass
        self.parent_point_pos = parent_point_pos

        if child_point_pos is None:
            child_point_pos = (0, 0, 0)  # Center of mass
        self.child_point_pos = child_point_pos

        self.parent_joint_vector = self._convert_tuple_to_vector(
            self.parent.get_frame(),
            self.parent_point_pos)
        self.child_joint_vector = self._convert_tuple_to_vector(
            self.child.get_frame(),
            self.child_point_pos)

        self._set_parent_child_rel()
        self._locate_joint_point()
        self.apply_joint()

    def _set_parent_child_rel(self):
        self.child.parent = self.parent
        self.parent.child = self.child

    def _convert_tuple_to_vector(self, frame, pos_tuple):
        if len(pos_tuple) != 3:
            raise TypeError('position tuple must be of length 3')
        else:
            unit_vectors = [frame.x, frame.y, frame.z]
            vector = Vector(0)
            for i in range(3):
                vector += pos_tuple[i] * unit_vectors[i]
            return vector

    def _locate_joint_point(self):
        self.parent_joint_point = self.parent.get_masscenter().locatenew(
            self.name + '_parent_joint',
            self.parent_joint_vector)
        self.child_joint_point = self.child.get_masscenter().locatenew(
            self.name + '_child_joint',
            self.child_joint_vector)

    def get_angle(self, vec1, vec2):
        mag1 = vec1.magnitude()
        mag2 = vec2.magnitude()
        return acos(dot(vec1, vec2)/(mag1 * mag2))

    def add_kd(self, kd):
        self._kd.append(kd)

    def get_kds(self):
        return self._kd

    def add_coordinate(self, coordinate):
        self._coordinates.append(coordinate)

    def get_coordinates(self):
        return self._coordinates

    def add_speed(self, speed):
        self._speed.append(speed)

    def get_speeds(self):
        return self._speed

    def apply_joint(self):
        """To create a custom joint, define a subclass of Joint class and
        and override this method"""
        raise NotImplementedError("To define a custom pydy.Joint, you need to" +
                                  " override apply_joint method in Joint's" +
                                  " subclass.")


class PinJoint(Joint):
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None, parent_axis=None, child_axis=None):

        if parent_axis is None:
            self.parent_axis = parent.get_frame().x
        else:
            self.parent_axis = parent_axis

        if child_axis is None:
            self.child_axis = child.get_frame().x
        else:
            self.child_axis = child_axis

        super(PinJoint, self).__init__(name, parent, child, parent_point_pos,
                                       child_point_pos)

    def align_axes(self):
        """Rotates child_frame such that child_axis is aligned to parent_axis.
        """
        angle = self.get_angle(self.parent_axis, self.child_axis)
        axis = cross(self.child_axis, self.parent_axis)
        if axis != Vector(0):
            self.child.get_frame().orient(
                self.parent.get_frame(), 'Axis', [angle, axis])

    def apply_joint(self):
        theta = dynamicsymbols('theta')
        thetad = dynamicsymbols('theta', 1)
        omega = dynamicsymbols('omega')
        self.add_coordinate(theta)
        self.add_speed(omega)
        self.add_kd(thetad - omega)
        self.child.get_frame().orient(self.parent.get_frame(), 'Axis',
                                      [theta, self.parent_axis])
        self.child.get_frame().set_ang_vel(self.parent.get_frame(), omega * self.parent_axis)
        self.align_axes()
        self.child_joint_point.set_pos(self.parent_joint_point, 0)
        self.child.get_masscenter().v2pt_theory(self.parent.get_masscenter(),
                                                self.parent.get_frame(),
                                                self.child.get_frame())


class SlidingJoint(Joint):
    def __init__(self, name, parent, child, parent_point_pos, child_point_pos,
                 parent_dir=None, child_dir=None):

        if parent_dir is None:
            self.parent_direction = parent.get_frame().x
        else:
            self.parent_direction = parent_dir

        if child_dir is None:
            self.child_direction = child.get_frame().x
        else:
            self.child_direction = child_dir

        super(SlidingJoint, self).__init__(name, parent, child, parent_point_pos,
                                           child_point_pos)

    def align_directions(self):
        """Rotates child_frame such that child_axis is aligned to parent_axis.
        """
        angle = self.get_angle(self.parent_direction, self.child_direction)
        if angle is not sympify(0):
            axis = cross(self.child_direction, self.parent_direction)
            self.child.get_frame().orient(
                self.parent.get_frame(), 'Axis', [angle, axis])

    def apply_joint(self):
        dis = dynamicsymbols('dis')
        disd = dynamicsymbols('dis', 1)
        vel = dynamicsymbols('vel')
        self.add_coordinate(dis)
        self.add_speed(vel)
        self.add_kd([disd - vel])
        self.child.get_frame().orient(self.parent.get_frame(), 'Axis',
                                      [0, self.parent.get_frame().z])
        self.align_directions()
        self._locate_joint_point()

        self.parent_joint_point.set_vel(self.parent.get_frame(), 0)
        self.child_joint_point.set_vel(self.child.get_frame(), 0)

        self.child_joint_point.set_pos(self.parent_joint_point,
                                       dis * self.parent_direction)
        self.child_joint_point.set_vel(self.parent.get_frame(),
                                       vel * self.parent_direction)
        self.child.get_masscenter().set_vel(self.parent.get_frame(),
                                            vel * self.parent_direction)
