from sympy.physics.mechanics import KanesMethod

__all__ = ['JointsMethod']


class JointsMethod(object):
    """
    # TODO
    Parameters
    ----------
    joints: list of Joints
        All the joints created must be passed as argument.
    root_body: Body
        Root body w.r.t which equations of motion much be generated.
        Is necessary to remove the  ambiguity since user may not
        connect all bodies using joints and thus, there is not
        way to find the root_body.
    """
    def __init__(self, joints, root_body):
        self.joints = joints
        self.root_body = root_body

    def get_all_bodies(self):
        self.bodies = []
        for joint in self.joints:
            parent = joint.parent
            if parent not in self.bodies:
                self.bodies.append(parent)
            child = joint.child
            if child not in self.bodies:
                self.bodies.append(child)

    def get_force_list(self):
        self.force_list = []
        for body in self.bodies:
            for force in body.force_list:
                self.force_list.append(force)

    def get_joints_details(self):
        self.kd = []
        self.q_ind = []
        self.u_ind = []
        for joint in self.joints:
            self.kd.append(joint.kd)
            self.q_ind.append(joint.q_ind)
            self.u_ind.append(joint.u_ind)

    def get_kanes(self):
        self.get_all_bodies()
        self.get_force_list()
        self.get_joints_details()
        self.KM = KanesMethod(self.root_body, q_ind=self.q_ind, u_ind=self.u_ind,
                              kd_eqs=self.kd)
        self.KM.kanes_equations(self.force_list, self.bodies)
        return self.KM
