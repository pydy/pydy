#!/usr/bin/env python

from sympy import Dummy, lambdify
from sympy.external import import_module

np = import_module('numpy')

def numeric_right_hand_side(kane, parameters):
    """Returns the right hand side of the first order ordinary differential
    equations from a KanesMethod system which can be evaluated numerically.

    This function probably only works for simple cases (no dependent speeds,
    specified functions).

    """

    dynamic = kane._q + kane._u
    dummy_symbols = [Dummy() for i in dynamic]
    dummy_dict = dict(zip(dynamic, dummy_symbols))
    kindiff_dict = kane.kindiffdict()

    M = kane.mass_matrix_full.subs(kindiff_dict).subs(dummy_dict)
    F = kane.forcing_full.subs(kindiff_dict).subs(dummy_dict)

    M_func = lambdify(dummy_symbols + parameters, M)
    F_func = lambdify(dummy_symbols + parameters, F)

    def right_hand_side(x, t, args):
        """Returns the derivatives of the states.

        Parameters
        ----------
        x : ndarray, shape(n)
            The current state vector.
        t : float
            The current time.
        args : ndarray
            The constants.

        Returns
        -------
        dx : ndarray, shape(2 * (n + 1))
            The derivative of the state.

        """
        arguments = np.hstack((x, args))
        dx = np.array(np.linalg.solve(M_func(*arguments),
            F_func(*arguments))).T[0]

        return dx

    return right_hand_side