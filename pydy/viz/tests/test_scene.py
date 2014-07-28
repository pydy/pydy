#!/usr/bin/env python

import os

from numpy import array, linspace
from scipy.integrate import odeint
from sympy import symbols
import sympy.physics.mechanics as me

from pydy.codegen.code import generate_ode_function
from pydy.viz import Sphere, VisualizationFrame, Scene


def test_create_static_html():
    # derive simple system
    mass, stiffness, damping, gravity = symbols('m, k, c, g')
    position, speed = me.dynamicsymbols('x v')
    positiond = me.dynamicsymbols('x', 1)
    ceiling = me.ReferenceFrame('N')
    origin = me.Point('origin')
    origin.set_vel(ceiling, 0)
    center = origin.locatenew('center', position * ceiling.x)
    center.set_vel(ceiling, speed * ceiling.x)
    block = me.Particle('block', center, mass)
    kinematic_equations = [speed - positiond]
    total_force = mass * gravity - stiffness * position - damping * speed
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
    system = (mass_matrix, forcing_vector, constants, coordinates, speeds)

    # integrate eoms
    evaluate_ode = generate_ode_function(*system)
    x0 = array([0.1, -1.0])
    args = {'constants': array([1.0, 1.0, 0.2, 9.8])}
    t = linspace(0.0, 10.0, 100)
    y = odeint(evaluate_ode, x0, t, args=(args,))

    # create visualization
    sphere = Sphere()
    viz_frame = VisualizationFrame(ceiling, block, sphere)
    scene = Scene(ceiling, origin, viz_frame)
    scene.generate_visualization_json([position, speed], list(constants), y,
                                      args['constants'])

    # test static dir creation
    scene.create_static_html(overwrite=True)
    assert os.path.exists('static')
    assert os.path.exists('static/index.html')
    assert os.path.exists('static/data.json')

    # test static dir deletion
    scene.remove_static_html(force=True)
    assert not os.path.exists('static')
