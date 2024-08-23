#!/usr/bin/env python

import os
import shutil
import glob

import pytest

import numpy as np
from numpy.testing import assert_allclose
from sympy import symbols
import sympy.physics.mechanics as me

from ...system import System
from ..shapes import Sphere
from ..visualization_frame import VisualizationFrame
from ..camera import PerspectiveCamera, OrthoGraphicCamera
from ..light import PointLight
from ..scene import Scene
from ...utils import sympy_newer_than


class TestScene(object):

    def setup_method(self):
        """Setups a simple 1 DoF mass spring damper system visualization."""

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

        if sympy_newer_than('1.0'):
            kane.kanes_equations(particles, forces)
        else:
            kane.kanes_equations(forces, particles)


        self.sys = System(kane)
        self.sys.initial_conditions = {position: 0.1, speed: -1.0}
        self.sys.constants = {mass: 1.0, stiffness: 2.0, damping: 3.0,
                              gravity: 9.8}
        self.sys.times = np.linspace(0.0, 0.01, 2)

        sphere = Sphere()

        self.ref_frame = ceiling
        self.origin = origin
        self.viz_frame = VisualizationFrame(ceiling, block, sphere)
        self.viz_frame_sym_shape = VisualizationFrame(ceiling, block,
                                                      Sphere(radius=mass /
                                                             10.0))

    def test_init(self):

        # test minimal args
        scene = Scene(self.ref_frame, self.origin)

        assert scene.reference_frame == self.ref_frame
        assert scene.origin == self.origin
        assert scene.name == 'unnamed'
        # NOTE : There isn't any way to compare the actual camera and light
        # objects because they will be created with different instances of a
        # point. So this is the best bet for now:
        assert scene.cameras[0].name == 'DefaultCamera'
        assert scene.lights[0].name == 'DefaultLight'

        scene.visualization_frames = [self.viz_frame]
        assert scene.visualization_frames == [self.viz_frame]
        assert scene.system is None
        assert scene.times is None
        assert scene.constants is None
        assert scene.states_symbols is None
        assert scene.states_trajectories is None
        assert scene.frames_per_second == 30

        # test maximal args/kwargs
        custom_camera = PerspectiveCamera('my_camera', self.ref_frame,
                                          self.origin.locatenew(
                                              'cam_point', 20 *
                                              self.ref_frame.z))

        custom_light = PointLight('my_light', self.ref_frame,
                                  self.origin.locatenew('light_point', 20 *
                                                        self.ref_frame.y))

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      name='my_scene', cameras=[custom_camera],
                      lights=[custom_light], system=self.sys)

        assert scene.visualization_frames == [self.viz_frame]
        assert scene.name == 'my_scene'
        assert scene.cameras == [custom_camera]
        assert scene.lights == [custom_light]
        assert scene.system is self.sys

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      name='my_scene', cameras=[custom_camera],
                      lights=[custom_light], times=self.sys.times,
                      constants=self.sys.constants,
                      states_symbols=self.sys.states,
                      states_trajectories=self.sys.integrate())

        assert scene.system is None
        assert_allclose(scene.times, self.sys.times)
        assert scene.constants == self.sys.constants
        assert scene.states_symbols == self.sys.states
        assert scene.states_trajectories.shape == (2, 2)

    def test_setting_incompatible_attrs(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      times=self.sys.times)

        with pytest.raises(ValueError):
            scene.system = self.sys

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      system=self.sys)

        with pytest.raises(ValueError):
            scene.times = self.sys.times

    def test_clear_trajectories(self):
        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      times=self.sys.times, constants=self.sys.constants,
                      states_symbols=self.sys.states,
                      states_trajectories=self.sys.integrate())
        scene.clear_trajectories()
        assert scene.system is None
        assert scene.times is None
        assert scene.constants is None
        assert scene.states_symbols is None
        assert scene.states_trajectories is None

    def test_generate_simulation_dict(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      system=self.sys)

        light_id = id(scene.lights[0])
        camera_id = id(scene.cameras[0])
        viz_frame_id = id(scene.visualization_frames[0])

        scene._generate_simulation_dict()

        expected_dict = {viz_frame_id: [[1.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,
                                         0.0, 0.0, 1.0, 0.0,
                                         0.1, 0.0, 0.0, 1.0],
                                        [1.0, 0.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0, 0.0,
                                         0.0, 0.0, 1.0, 0.0,
                                         0.09062405469543587, 0.0, 0.0, 1.0]],
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

        assert (sorted(expected_dict.keys()) ==
                sorted(scene._simulation_info.keys()))

        for k, v in scene._simulation_info.items():
            assert_allclose(v, expected_dict[k])

    def test_generate_scene_dict(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      system=self.sys)

        light_id = id(scene.lights[0])
        camera_id = id(scene.cameras[0])
        viz_frame_id = id(scene.visualization_frames[0])

        # NOTE : generate_simulation_dict must be called before
        # generate_scene_dict
        scene._generate_simulation_dict()
        scene._generate_scene_dict()

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
                                     'color': 'grey',
                                     'material': 'default',
                                     'reference_frame_name': 'N',
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

        assert scene._scene_info == expected_scene_dict

        # Test the constant_map kwarg.
        scene = Scene(self.ref_frame, self.origin, self.viz_frame_sym_shape,
                      system=self.sys)
        viz_frame_id = id(scene.visualization_frames[0])
        scene._generate_simulation_dict()
        scene._generate_scene_dict()

        assert scene._scene_info['objects'][viz_frame_id]['radius'] == 0.1

    def test_custom_camera(self):

        camera_frame = self.ref_frame.orientnew('rot', 'body',
                                                [np.pi / 2.0,
                                                 np.pi / 2.0,
                                                 np.pi / 2.0], 'xyz')

        camera_point = self.origin.locatenew('cam_point', 30.0 *
                                             camera_frame.z)

        camera = OrthoGraphicCamera('my_camera', camera_frame, camera_point)

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      cameras=[camera], system=self.sys)

        camera_id = id(camera)

        scene._generate_simulation_dict()
        scene._generate_scene_dict()
        scene_dict = scene._scene_info

        expected_orientation_matrix = np.array([0.0, 0.0, 1.0, 0.0,
                                                0.0, -1.0, 0.0, 0.0,
                                                1.0, 0.0, 0.0, 0.0,
                                                30.0, 0.0, 0.0, 1.0])

        assert_allclose(scene_dict['cameras'][camera_id]['init_orientation'],
                        expected_orientation_matrix, atol=1e-14)

        assert scene_dict['cameras'][camera_id]['type'] == 'OrthoGraphicCamera'

    def test_create_static_html(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame,
                      system=self.sys)

        # test static dir creation
        scene.create_static_html(overwrite=True, silent=True, prefix='test')
        assert os.path.exists('pydy-resources')
        assert os.path.exists('pydy-resources/index.html')
        assert os.path.exists('pydy-resources/test_scene_desc.json')
        assert os.path.exists('pydy-resources/test_simulation_data.json')

        # test static dir deletion
        scene.remove_static_html(force=True)
        assert not os.path.exists(Scene.pydy_directory)

    def test_generate_visualization_json_system(self):

        scene = Scene(self.ref_frame, self.origin, self.viz_frame)
        scene.generate_visualization_json_system(self.sys)

        # Tests issue #204
        assert scene._scene_info['constant_map'] == {'m': 1.0, 'k': 2.0,
                                                     'c': 3.0, 'g': 9.8}

    def teardown_method(self):

        try:
            shutil.rmtree(Scene.pydy_directory)
        except OSError:
            pass

        for json in glob.glob("*.json"):
            os.remove(json)
