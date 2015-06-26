__all__ = ['PinJoint', 'RevoluteJoint']

class Joint(object):
    def __init__(self, parent, child):
        self.parent = parent
        self.child = child
        # TODO

class PinJoint(Joint):
    def __init__(self, parent, parent_vector, child, child_vector, axis=None):
        Joint.__init__(self, parent, child)
        # TODO

class RevoluteJoint(Joint):
    def __init__(self, parent, child, **kwargs):
        Joint.__init__(self, parent, child)
        # TODO
