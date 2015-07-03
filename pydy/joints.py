__all__ = ['PinJoint', 'PrismaticJoint', 'CylindricJoint', 'SphericalJoint', 'PlanarJoint']


class Joint(object):
    """Serves as a base class for all specific Joints. Contains all public methods
    to be used by other specific Joint's implementation.

    Note: Some of these methods can be shifted to body. Want to have a discussion
    on this.
    """
    def __init__(self, name, parent, parent_joint_vector, child, child_joint_vector):
        self.name = name
        self.parent = parent
        self.child = child
        self.parent_joint_vector = parent_joint_vector
        self.child_joint_vector = child_joint_vector

    def set_parent_child_rel(self):
        self.child.parent = self.parent
        self.parent.child = self.child

    def _locate_joint_point(self):
        """Locates a new point, Point of Joint, using parent's center of mass
        and parent_joint_vector and sets position of child using child_joint_vector
        from joint point."""
        self.joint_point = self.parent.masscenter.locatenew(self.name + '_Point', self.parent_joint_vector)
        self.child.masscenter.set_pos(self.joint_point, self.child_joint_vector)

    def set_joint_point_vel(self, value):
        """Sets the velocity of the point at joint w.r.t parent's frame."""
        # TODO

    def set_child_vel(self, value):
        """sets the velocity of child masscenter to value w.r.t parent's frame.
        Note that the velocity is added to point where as angular velocity is
        added to frame."""
        # TODO

    def set_child_angular_vel(self, value):
        """sets angular velocity of child's frame w.r.t. to parent's frame.
        Note that the velocity is added to point where as angular velocity is
        added to frame."""
        # TODO

    def orient_child(self):
        # Note: Can be removed and just orient() be used, including it
        # here as a reference of available functionalities.
        """Orients child's frame w.r.t parent's frame. Please refer to
        sympy.physics.vector.frame.ReferenceFrame.orient() for all the
        paramter details. Its exactly same"""
        # TODO

    def _apply_joint(self):
        """For custom joint, make a subclass of Joint and override this method."""
        #TODO generic joint creation.


class PinJoint(Joint):
    """Uses methods in Joint's class and create a Revolute (Pin) Joint between
    parent and child."""
    def __init__(self, name, parent, parent_joint_vector, child, child_joint_vector, axis=None):
        super(Joint, self).__init__(*args, **kwargs)
        if axis is None:
            self.axis = parent.frame.z
        else:
            self.axis = axis
        self._apply_joint()

    def _apply_joint(self):
        """Specific implementation of a joint using all public and private methods."""
        #TODO


class PrismaticJoint(Joint):
    def __init__(self, name, parent, parent_joint_vector, child, child_joint_vector, direction=None):
        super(Joint, self).__init__(*args, **kwargs)
        if direction is None:
            self.direction = self.parent.z
        else:
            self.direction = direction
        self._apply_joint()

    def _apply_joint(self):
        """Specific implementation of a joint using all public and private methods."""
        #TODO

# Other classes will have similar implementation

class CylindricJoint(Joint):
    #TODO

class SphericalJoint(Joint):
    #TODO

class PlanarJoint(Joint):
    #TODO
