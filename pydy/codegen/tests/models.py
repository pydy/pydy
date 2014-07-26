#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains some sample symbolic models used for testing and
examples."""

# external libraries
from sympy import symbols
import sympy.physics.mechanics as me


def generate_mass_spring_damper_equations_of_motion(external_force=True):
    """Returns the symbolic equations of motion and associated variables for a
    simple one degree of freedom mass, spring, damper system with gravity
    and an optional external specified force.

      / / / / / / / /
     ----------------
      |    |     |   | g
      \   | |    |   V
    k /   --- c  |
      |    |     | x, v
     --------    V
     |  m   | -----
     --------
        | F
        V

    Returns
    -------
    mass_matrix : sympy.MutableMatrix, shape(2,2)
        The symbolic mass matrix of the system which are linear in
        derivatives of the states.
    forcing_vector : sympy.MutableMatrix, shape(2,1)
        The forcing vector of the system.
    constants : tuple, len(4)
        A sequence of all the symbols which are constants in the equations
        of motion, i.e. m, k, c, g.
    coordinates : tuple, len(1)
        A sequence of all the dynamic symbols, i.e. functions of time, which
        describe the configuration of the system, i.e. x.
    speeds : tuple, len(1)
        A sequence of all the dynamic symbols, i.e. functions of time, which
        describe the generalized speeds of the system, i.e v.
    specified : tuple, len(1) or None
        A sequence of all the dynamic symbols, i.e. functions of time, which
        describe the specified variables of the system, i.e F. If
        `external_force` is False, then None is returned.
    external_force : boolean, default=True
        If true a specifed force will be added to the system.

    """

    mass, stiffness, damping, gravity = symbols('m, k, c, g')

    position, speed = me.dynamicsymbols('x v')
    positiond = me.dynamicsymbols('x', 1)
    force = me.dynamicsymbols('F')

    ceiling = me.ReferenceFrame('N')

    origin = me.Point('origin')
    origin.set_vel(ceiling, 0)

    center = origin.locatenew('center', position * ceiling.x)
    center.set_vel(ceiling, speed * ceiling.x)

    block = me.Particle('block', center, mass)

    kinematic_equations = [speed - positiond]

    total_force = mass * gravity - stiffness * position - damping * speed
    if external_force is True:
        total_force += force
    forces = [(center, total_force * ceiling.x)]

    particles = [block]

    kane = me.KanesMethod(ceiling, q_ind=[position], u_ind=[speed],
                          kd_eqs=kinematic_equations)
    kane.kanes_equations(forces, particles)

    mass_matrix = kane.mass_matrix_full
    forcing_vector = kane.forcing_full
    constants = (mass, stiffness, damping, gravity)
    coordinates = (position,)
    speeds = (speed,)
    if external_force is True:
        specified = (force,)
    else:
        specified = None

    return (mass_matrix, forcing_vector, constants, coordinates, speeds,
            specified)


def generate_n_link_pendulum_on_cart_equations_of_motion(n, cart_force=True,
                                                         joint_torques=False):
    """Returns the the symbolic first order equations of motion for a 2D
    n-link pendulum on a sliding cart under the influence of gravity in this
    form:

        M(x) x(t) = F(x, u, t)

    Parameters
    ----------
    n : integer
        The number of links in the pendulum.
    cart_force : boolean, default=True
        If true an external specified lateral force is applied to the cart.
    joint_torques : boolean, default=False
        If true joint torques will be added as specified inputs at each
        joint.

    Returns
    -------
    mass_matrix : sympy.MutableMatrix, shape(2 * (n + 1), 2 * (n + 1))
        The symbolic mass matrix of the system which are linear in u' and q'.
    forcing_vector : sympy.MutableMatrix, shape(2 * (n + 1), 1)
        The forcing vector of the system.
    constants : list
        A sequence of all the symbols which are constants in the equations
        of motion.
    coordinates : list
        A sequence of all the dynamic symbols, i.e. functions of time, which
        describe the configuration of the system.
    speeds : list
        A sequence of all the dynamic symbols, i.e. functions of time, which
        describe the generalized speeds of the system.
    specfied : list
        A sequence of all the dynamic symbols, i.e. functions of time, which
        describe the specified inputs to the system.

    Notes
    -----
    The degrees of freedom of the system are n + 1, i.e. one for each
    pendulum link and one for the lateral motion of the cart.

    M x' = F, where x = [u0, ..., un+1, q0, ..., qn+1]

    The joint angles are all defined relative to the ground where the x axis
    defines the ground line and the y axis points up. The joint torques are
    applied between each adjacent link and the between the cart and the
    lower link where a positive torque corresponds to positive angle.

    """
    if n <= 0:
        raise ValueError('The number of links must be a positive integer.')

    q = me.dynamicsymbols('q:{}'.format(n + 1))
    u = me.dynamicsymbols('u:{}'.format(n + 1))

    if joint_torques is True:
        T = me.dynamicsymbols('T1:{}'.format(n + 1))

    m = symbols('m:{}'.format(n + 1))
    l = symbols('l:{}'.format(n))
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

    if cart_force is True or joint_torques is True:
        specified = []
    else:
        specified = None

    for i in range(n):
        Bi = I.orientnew('B{}'.format(i), 'Axis', [q[i + 1], I.z])
        Bi.set_ang_vel(I, u[i + 1] * I.z)
        frames.append(Bi)

        Pi = points[-1].locatenew('P{}'.format(i + 1), l[i] * Bi.x)
        Pi.v2pt_theory(points[-1], I, Bi)
        points.append(Pi)

        Pai = me.Particle('Pa' + str(i + 1), Pi, m[i + 1])
        particles.append(Pai)

        forces.append((Pi, -m[i + 1] * g * I.y))

        if joint_torques is True:

            specified.append(T[i])

            if i == 0:
                forces.append((I, -T[i] * I.z))

            if i == n - 1:
                forces.append((Bi, T[i] * I.z))
            else:
                forces.append((Bi, T[i] * I.z - T[i + 1] * I.z))

        kindiffs.append(q[i + 1].diff(t) - u[i + 1])

    if cart_force is True:
        F = me.dynamicsymbols('F')
        forces.append((P0, F * I.x))
        specified.append(F)

    kane = me.KanesMethod(I, q_ind=q, u_ind=u, kd_eqs=kindiffs)
    kane.kanes_equations(forces, particles)

    mass_matrix = kane.mass_matrix_full
    forcing_vector = kane.forcing_full
    coordinates = kane._q
    speeds = kane._u

    constants = [g, m[0]]
    for i in range(n):
        constants += [l[i], m[i + 1]]

    return (mass_matrix, forcing_vector, constants, coordinates, speeds,
            specified)
