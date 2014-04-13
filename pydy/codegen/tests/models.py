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


def generate_n_link_pendulum_on_cart_equations_of_motion(n, cart_force=True):
    """Returns the implicit form of the symbolic equations of motion for a
    2D n-link pendulum on a sliding cart under the influence of gravity.

    Parameters
    ----------
    n : integer
        The number of links in the pendulum.
    cart_force : boolean, default=True
        If true an external specified lateral force is applied to the cart.

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


    Notes
    -----
    The degrees of freedom of the system are n + 1, i.e. one for each
    pendulum link and one for the lateral motion of the cart.

    M x' = F, where x = [u0, ..., un+1, q0, ..., qn+1]

    """

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

    if cart_force is True:
        F = me.dynamicsymbols('F')
        forces.append((P0, F * I.x))
        specified = [F]
    else:
        specified = None

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
