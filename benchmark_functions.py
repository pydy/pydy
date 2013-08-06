#!/usr/bin/env python

import time
from sympy import symbols, Dummy, lambdify
from sympy.utilities.autowrap import autowrap
import sympy.physics.mechanics as me
from sympy.printing.theanocode import theano_function
import numpy as np

# debugging
try:
    from IPython.core.debugger import Tracer
except ImportError:
    pass
else:
    set_trace = Tracer()

#class System():
    #track parameters and their values
    # u's and q's

def generate_equations(n):
    """Returns a KanesMethod object which contains the symbolic equations of
    motion for a 2D n-pendulum on a sliding cart."""

    q = me.dynamicsymbols('q:' + str(n + 1))
    u = me.dynamicsymbols('u:' + str(n + 1))

    m = symbols('m:' + str(n + 1))
    l = symbols('l:' + str(n))
    g, t = symbols('g t')

    I = me.ReferenceFrame('I')
    O = me.Point('O')
    O.set_vel(I, 0)

    P0 = me.Point('P0')
    P0.set_pos(O, q[0] * I.x)
    P0.set_vel(I, u[0] * I.x)
    Pa0 = me.Particle('Pa0', P0, m[0])

    frames = [I]
    points = [P0]
    particles = [Pa0]
    forces = [(P0, -m[0] * g * I.y)]
    kindiffs = [q[0].diff(t) - u[0]]

    for i in range(n):
        Bi = I.orientnew('B' + str(i), 'Axis', [q[i + 1], I.z])
        Bi.set_ang_vel(I, u[i + 1] * I.z)
        frames.append(Bi)

        Pi = points[-1].locatenew('P' + str(i + 1), l[i] * Bi.x)
        Pi.v2pt_theory(points[-1], I, Bi)
        points.append(Pi)

        Pai = me.Particle('Pa' + str(i + 1), Pi, m[i + 1])
        particles.append(Pai)

        forces.append((Pi, -m[i + 1] * g * I.y))

        kindiffs.append(q[i + 1].diff(t) - u[i + 1])

    kane = me.KanesMethod(I, q_ind=q, u_ind=u, kd_eqs=kindiffs)
    kane.kanes_equations(forces, particles)

    parameters = [g, m[0]]
    for i in range(n):
        parameters += [l[i], m[i + 1]]

    return kane, parameters

def generate_numeric_eom_matrices(mass_matrix, forcing_vector, constants,
        states, specified=None, generator="lambdify", **options):
    """Returns functions which compute the mass matrix and forcing vector given
    numerical values of the states, constants, and optionally specified
    variables.

    Parameters
    ----------
    mass_matrix : SymPy Matrix, shape(n, n)
        An n x n symbolic matrix where each entry is an expression made up of
        the states, constants, and specified variables.
    forcing_vector : SymPy Matrix, shape(n,)
        An n x 1 symbolic matrix where each entry is an expression made up of
        the states, constants, and specified variables.
    constants : sequence of sympy.Symbol
    states : sequence of sympy.Function
    specified : sequence of sympy.Function
    generator : string, optional, default='lambdify'

    Returns
    -------
    mass_matrix_func : function
    forcing_vector_func : function

    constants + states [+ specified]

    Examples
    --------

    >>> q = dynamicsymbols('q:2')
    >>> u = dynamicsymbols('u:2')
    >>> f = dynamicsymbols('f:2')
    >>> c = symbols('c:2')
    >>> mass_matrix = Matrix([[q[0] + u[0] + c[0], c[0] * (q[1] + u[1])],

    >>> type(mass_matrix)
    sympy.matrices.dense.MutableDenseMatrix
    >>> type(forcing_vector)
    sympy.matrices.dense.MutableDenseMatrix
    >>> generate_numeric_matrices(mass_matrix, forcing_vector, c, q + u,
        >>> specified=f,


    """

    dynamic = states
    if specified is not None:
        dynamic += specified

    # TODO : the Dummy symbols may not be needed now that lambdify is updated.
    # the theano option may not need them either.
    dummy_symbols = [Dummy() for i in dynamic]
    dummy_dict = dict(zip(dynamic, dummy_symbols))

    dummy_mass_matrix = mass_matrix.subs(dummy_dict)
    dummy_forcing_vector = forcing_vector.subs(dummy_dict)

    arguments = constants + dummy_symbols

    if generator == 'lambdify':
        mass_matrix_func = lambdify(arguments, dummy_mass_matrix)
        forcing_vector_func = lambdify(arguments, dummy_forcing_vector)
    elif generator == 'theano':
        mass_matrix_func = theano_function(arguments, [dummy_mass_matrix],
                                           on_unused_input='ignore')
        forcing_vector_func = theano_function(arguments,
                                              [dummy_forcing_vector],
                                              on_unused_input='ignore')
        # lower run time from 0.15s to 0.08s for n=1
        mass_matrix_func.trust_input = True
        forcing_vector_func.trust_input = True

    elif generator == 'autowrap':
        funcs = []
        for entry in dummy_mass_matrix:
            funcs.append(autowrap(entry, args=arguments, **options))

        def mass_matrix_func(*args):
            result = []
            for func in funcs:
                result.append(func(*args))
            # TODO : this may not be correctly reshaped
            return np.matrix(result).reshape(np.sqrt(len(result)),
                                             np.sqrt(len(result)))

        funcs = []
        for row in dummy_forcing_vector:
            funcs.append(autowrap(row, args=arguments, **options))

        def forcing_vector_func(*args):
            result = []
            for func in funcs:
                result.append(func(*args))
            return np.matrix(result)
    else:
        raise NotImplementedError('{} is not implemented yet'.format(generator))

    return mass_matrix_func, forcing_vector_func


def numeric_right_hand_side(kane, parameters, specified=None, generator='lambdify'):
    """Returns the right hand side of the first order ordinary differential
    equations from a KanesMethod system which can be evaluated numerically.

    This function probably only works for simple cases (no dependent speeds,
    specified functions, etc).

    Parameters
    ----------
    kane : KanesMethod object
    parameters : a sequence
        A list of symbols for the system constants.
    specified : a sequence
        A sequence of symbols/functions for specifed variables.
    generator : optional
        lambdify theano numba cython fortran parakeet

    Returns
    -------

    """

    dynamic = kane._q + kane._u

    kindiff_dict = kane.kindiffdict()

    mass_matrix = kane.mass_matrix_full.subs(kindiff_dict)
    forcing_vector = kane.forcing_full.subs(kindiff_dict)

    mass_matrix_func, forcing_vector_func = \
        generate_numeric_eom_matrices(mass_matrix, forcing_vector,
                                      parameters, dynamic,
                                      specified=specified,
                                      generator=generator)

    arguments = parameters + dynamic
    if specified is not None: arguments += specified
    start = time.time()
    numtimes = 1000

    # Lower from 0.5s to 0.15s for the run time when n=1
    #inp = np.random.random(len(arguments))
    inp = [np.asarray(x) for x in np.random.random(len(arguments))]
    for i in range(numtimes):
        mass_matrix_func(*inp)
        forcing_vector_func(*inp)
    total = time.time() - start
    print("It took {} seconds to compute M and F with {} {} times at an \
average of {} seconds per computation.".format(total, generator, numtimes,
    total / numtimes))

    #set_trace()

    def right_hand_side(x, t, args):
        """Returns the derivatives of the states.

        Parameters
        ----------
        x : ndarray, shape(n)
            The current state vector.
        t : float
            The current time.
        args : ndarray
            The specified variables and constants.

        Returns
        -------
        dx : ndarray, shape(2 * (n + 1))
            The derivative of the state.

        """
        arguments = np.hstack((x, args))
        dx = np.array(np.linalg.solve(mass_matrix_func(*arguments),
                                      forcing_vector_func(*arguments))).T[0]

        return dx

    return right_hand_side
