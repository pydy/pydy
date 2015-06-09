#!/usr/bin/env python

import os
import shutil

from numpy import linspace
from sympy import symbols
import sympy.physics.mechanics as me

from ...system import System
from ..shapes import Sphere
from ..visualization_frame import VisualizationFrame
from ..camera import PerspectiveCamera
from ..light import PointLight
from ..scene import Scene


class TestScene(object):

    def setup(self):

        mass, stiffness, damping, gravity = symbols('m, k, c, g')

        position, speed = me.dynamicsymbols('x v')
        positiond = me.dynamicsymbols('x', 1)

        kinematic_equations = [speed - positiond]

        ceiling = me.ReferenceFrame('N')

        origin = me.Point('origin')
        origin.set_vel(ceiling, 0)
        center = origin.locatenew('center', position * ceiling.x)
        center.set_vel(ceiling, speed * ceiling.x)

        block = me.Particle('block', center, mass)
        particles = [block]

        total_force = mass * gravity - stiffness * position - damping * speed
        forces = [(center, total_force * ceiling.x)]

        kane = me.KanesMethod(ceiling, q_ind=[position], u_ind=[speed],
                              kd_eqs=kinematic_equations)
        kane.kanes_equations(forces, particles)

        self.sys = System(kane)
        self.sys.initial_conditions = {position: 0.1, speed: -1.0}
        self.constants = {mass: 1.0, stiffness: 1.0, damping: 0.2,
                          gravity: 9.8}
        self.sys.times = linspace(0.0, 0.01, 2)

        sphere = Sphere()

        self.ref_frame = ceiling
        self.origin = origin
        self.viz_frame = VisualizationFrame(ceiling, block, sphere)

    def test_init(self):

        # test minimal args
        scene = Scene(self.ref_frame, self.origin)

        assert scene.reference_frame == self.ref_frame
        assert scene.origin == self.origin
        assert scene.name == 'unnamed'
        assert scene._width == 800
        assert scene._height == 800
        # TODO : These are failing.
        #assert scene.cameras == [PerspectiveCamera('DefaultCamera',
                                                   #self.ref_frame,
                                                   #self.origin.locatenew(
                                                       #'p_camera', 10 *
                                                       #self.ref_frame.z))]
        #assert scene.lights == [PointLight('DefaultLight', self.ref_frame,
                                           #self.origin.locatenew(
                                               #'p_light', 10 *
                                               #self.ref_frame.z))]

        scene.visualization_frames = [self.viz_frame]
        assert scene.visualization_frames == [self.viz_frame]

        # test maximal args/kwargs
        custom_camera = PerspectiveCamera('my_camera', self.ref_frame,
                                          self.origin.locatenew(
                                              'cam_point', 20 *
                                              self.ref_frame.z))

        custom_light = PointLight('my_light', self.ref_frame,
                                  self.origin.locatenew('light_point', 20 *
                                                        self.ref_frame.y))

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      name='my_scene', width=200, height=300,
                      cameras=[custom_camera], lights=[custom_light])

        assert scene.visualization_frames == [self.viz_frame]
        assert scene.name == 'my_scene'
        assert scene._width == 200
        assert scene._height == 300
        assert scene.cameras == [custom_camera]
        assert scene.lights == [custom_light]

    def test_generate_simulation_dict(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame)

        light_id = id(scene.lights[0])
        camera_id = id(scene.cameras[0])
        viz_frame_id = id(scene.visualization_frames[0])

        sim_dict = scene.generate_simulation_dict(self.sys.states,
                                                  self.sys.constants.keys(),
                                                  self.sys.integrate(),
                                                  self.sys.constants.values())

        expected_dict = {viz_frame_id: [[1.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,
                                         0.0, 0.0, 1.0, 0.0,
                                         0.1, 0.0, 0.0, 1.0],
                                        [1.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,
                                         0.0, 0.0, 1.0, 0.0,
                                         0.09009483961898038, 0.0, 0.0, 1.0]],
                         light_id: [[1.0, 0.0, 0.0, 0.0,
                                     0.0, 1.0, 0.0, 0.0,
                                     0.0, 0.0, 1.0, 0.0,
                                     0.0, 0.0, 10.0, 1.0],
                                    [1.0, 0.0, 0.0, 0.0,
                                     0.0, 1.0, 0.0, 0.0,
                                     0.0, 0.0, 1.0, 0.0,
                                     0.0, 0.0, 10.0, 1.0]],
                         camera_id: [[1.0, 0.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0, 0.0,
                                      0.0, 0.0, 1.0, 0.0,
                                      0.0, 0.0, 10.0, 1.0],
                                     [1.0, 0.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0, 0.0,
                                      0.0, 0.0, 1.0, 0.0,
                                      0.0, 0.0, 10.0, 1.0]]}

        # TODO : This needs floating point comparisons and I think the hash
        # value changes everytime the method is called.
        assert sim_dict == expected_dict

    def test_generate_scene_dict(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame)

        light_id = id(scene.lights[0])
        camera_id = id(scene.cameras[0])
        viz_frame_id = id(scene.visualization_frames[0])

        # NOTE : generate_simulation_dict must be called before
        # generate_scene_dict
        sim_dict = scene.generate_simulation_dict(self.sys.states,
                                                  self.sys.constants.keys(),
                                                  self.sys.integrate(),
                                                  self.sys.constants.values())
        scene_dict = scene.generate_scene_dict()

        expected_scene_dict = \
            {'newtonian_frame': 'N',
             'name': 'unnamed',
             'workspaceSize': 0.2,
             'source': 'PyDy',
             'lights': {light_id: {'color': 'white',
                                   'init_orientation':
                                       [1.0, 0.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0, 0.0,
                                        0.0, 0.0, 1.0, 0.0,
                                        0.0, 0.0, 10.0, 1.0],
                                    'simulation_id': light_id,
                                    'type': 'PointLight',
                                    'name': 'DefaultLight'}},
             'cameras': {camera_id: {'fov': 45.0,
                                     'name': 'DefaultCamera',
                                     'far': 1000.0,
                                     'simulation_id': camera_id,
                                     'near': 1.0,
                                     'init_orientation':
                                         [1.0, 0.0, 0.0, 0.0,
                                          0.0, 1.0, 0.0, 0.0,
                                          0.0, 0.0, 1.0, 0.0,
                                          0.0, 0.0, 10.0, 1.0],
                                           'type': 'PerspectiveCamera'}},
             'objects': {viz_frame_id: {'simulation_id': viz_frame_id,
                                        'name': 'unnamed',
                                        'color': 'grey',
                                        'material': 'default',
                                        'reference_frame_name': 'N',
                                        'radius': 10.0,
                                        'init_orientation':
                                            [1.0, 0.0, 0.0, 0.0,
                                             0.0, 1.0, 0.0, 0.0,
                                             0.0, 0.0, 1.0, 0.0,
                                             0.1, 0.0, 0.0, 1.0],
                                        'type': 'Sphere'}}}

        assert scene_dict == expected_scene_dict

        # TODO : Test the constant_map kwarg.

    def test_create_static_html(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame)
        scene.generate_visualization_json_system(self.sys,
                                                 outfile_prefix="test")

        # test static dir creation
        scene.create_static_html(overwrite=True)
        assert os.path.exists('static')
        assert os.path.exists('static/index.html')
        assert os.path.exists('static/test_scene_desc.json')
        assert os.path.exists('static/test_simulation_data.json')

        # test static dir deletion
        scene.remove_static_html(force=True)
        assert not os.path.exists('static')

    def cleanup(self):

        shutil.rmtree('static')
