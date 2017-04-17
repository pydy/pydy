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
        self._set_kanes()

    @property
    def bodylist(self):
        bodies = []
        body = self.root_body
        while body:
            bodies.append(body)
            body = body.child
        return bodies

    @property
    def forcelist(self):
        force_list = []
        for body in self.bodylist:
            for force in body.force_list:
                force_list.append(force)
        return force_list

    @property
    def q(self):
        q_ind = []
        for joint in self.joints:
            coordinates = joint.get_coordinates()
            for coordinate in coordinates:
                q_ind.append(coordinate)
        return q_ind

    @property
    def u(self):
        u_ind = []
        for joint in self.joints:
            speeds = joint.get_speeds()
            for speed in speeds:
                u_ind.append(speed)
        return u_ind

    @property
    def kd(self):
        kd_ind = []
        for joint in self.joints:
            kds = joint.get_kds()
            for kd in kds:
                kd_ind.append(kd)
        return kd_ind

    @property
    def forcing_full(self):
        return self._KM.forcing_full

    @property
    def forcing(self):
        return self._KM.forcing

    @property
    def mass_matrix_full(self):
        return self._KM.mass_matrix_full

    @property
    def mass_matrix(self):
        return self._KM.mass_matrix

    @property
    def auxiliary_eqs(self):
        return self._KM.auxiliary_eqs

    def _set_kanes(self):
        self._KM = KanesMethod(self.root_body.get_frame(), q_ind=self.q, u_ind=self.u,
                               kd_eqs=self.kd)
        self._KM.kanes_equations(self.forcelist, self.bodylist)
        # TODO Removing call to private attributes in pydy.System and fix this.
        self._qdot = self._KM._qdot
        self._udot = self._KM._udot
        self._uaux = self._KM._uaux
