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

        if parent_point_pos is None:
            parent_point_pos = (0,0,0)  # Center of mass
        self.parent_point_pos = parent_point_pos

        if child_point_pos is None:
            child_point_pos = (0,0,0)  # Center of mass
        self.child_point_pos = child_point_pos

        self.parent_joint_vector = self._convert_tuple_to_vector(
            self.parent.frame,
            self.parent_point_pos)
        self.child_joint_vector = self._convert_tuple_to_vector(
            self.child.frame,
            self.child_point_pos)

        self._set_parent_child_rel()
        self.apply_joint()

    def _set_parent_child_rel(self):
        self.child.parent = self.parent
        self.parent.child = self.child

    def _convert_tuple_to_vector(self, frame, pos_tuple):
        if len(pos_tuple) == 3:
            unit_vectors = [frame.x, frame.y, frame.z]
            vector = 0
            for i in range(3):
                vector += pos_tuple[i] * unit_vectors[i]
            return vector
        else:
            raise TypeError('position tuple must be of length 3')

    def _locate_joint_point(self):
        self.parent_joint_point = self.parent.masscenter.locatenew(
            self.name + '_parent_joint',
            self.parent_joint_vector)
        self.child_joint_point = self.child.masscenter.locatenew(
            self.name + '_child_joint',
            self.child_joint_vector)

    def apply_joint(self):
        """To create a custom joint, define a subclass of Joint class and
        and override this method"""
        raise NotImplementedError("To define a custom pydy.Joint, you need to" +
                                  " override apply_joint method in Joint's" +
                                  " subclass.")


class PinJoint(Joint):
    def __init__(self, name, parent, child, parent_point_pos=None,
                 child_point_pos=None):
        super(Joint, self).__init__()

    def apply_joint(self):
        # TODO


class SlidingJoint(Joint):
    def __init__(self, name, parent, child, parent_point_pos, child_point_pos):
        super(Joint, self).__init__()

    def apply_joint(self):
        # TODO


class CylindricJoint(Joint):
    def __init__(self, name, parent, child, parent_point_pos, child_point_pos):
        super(Joint, self).__init__()

    def apply_joint(self):
        # TODO


class SphericalJoint(Joint):
    def __init__(self, name, parent, child, parent_point_pos, child_point_pos):
        super(Joint, self).__init__()

    def apply_joint(self):
        # TODO


class PlanarJoint(Joint):
    def __init__(self, name, parent, child, parent_point_pos, child_point_pos):
        super(Joint, self).__init__()

    def apply_joint(self):
        # TODO
