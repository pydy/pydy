__all__ = ['JointsMethod']

class JointsMethod(object):
    """
    # TODO
    Parameters
    ----------
    root_body: Body
        Root body in a system where bodies are connected using joints.
    """
    def __init__(self, root_body):
        self.root_body = root_body

    def get_all_bodies(self):
        self.bodies = []
        body = self.root_body
        while(body):
            self.bodies.append(body)
            body = body.child

    def get_equations(self):
        self.get_all_bodies()
        # TODO get all values from bodies
        # TODO generate equations using Kane's method.




