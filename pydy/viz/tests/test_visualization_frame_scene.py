#!/usr/bin/env python

# external
from sympy import sin, cos, symbols, pi
from sympy.physics.mechanics import (dynamicsymbols, ReferenceFrame, Point,
                                     RigidBody, Particle, inertia)
from numpy import radians
from numpy.testing import assert_allclose

# local
from ..shapes import Cylinder
from ..visualization_frame import VisualizationFrame
from ..camera import PerspectiveCamera, OrthoGraphicCamera
from ..scene import Scene
from ..light import PointLight


class TestVisualizationFrameScene(object):

    def __init__(self):

        self.coords = dynamicsymbols('q:3')
        self.speeds = dynamicsymbols('p:3')

        self.dynamic = list(self.coords) + list(self.speeds)
        self.states = ([radians(45) for x in self.coords] +
                       [radians(30) for x in self.speeds])

        self.inertial_ref_frame = ReferenceFrame('I')
        self.A = self.inertial_ref_frame.orientnew('A', 'space',
                                                   self.coords, 'XYZ')
        self.B = self.A.orientnew('B', 'space', self.speeds, 'XYZ')
        self.camera_ref_frame = \
            self.inertial_ref_frame.orientnew('C', 'space', (pi / 2, 0, 0),
                                              'XYZ')

        self.origin = Point('O')
        self.P1 = self.origin.locatenew('P1',
                                        10 * self.inertial_ref_frame.x +
                                        10 * self.inertial_ref_frame.y +
                                        10 * self.inertial_ref_frame.z)
        self.P2 = self.P1.locatenew('P2',
                                    10 * self.inertial_ref_frame.x +
                                    10 * self.inertial_ref_frame.y +
                                    10 * self.inertial_ref_frame.z)

        self.point_list1 = [[2, 3, 1], [4, 6, 2], [5, 3, 1], [5, 3, 6]]
        self.point_list2 = [[3, 1, 4], [3, 8, 2], [2, 1, 6], [2, 1, 1]]

        self.shape1 = Cylinder(1.0, 1.0)
        self.shape2 = Cylinder(1.0, 1.0, name='booger')

        self.Ixx, self.Iyy, self.Izz = symbols('Ixx, Iyy, Izz')
        self.mass = symbols('m')
        self.parameters = [self.Ixx, self.Iyy, self.Izz, self.mass]
        self.param_vals = [2.0, 3.0, 4.0, 5.0]
        self.inertia = inertia(self.A, self.Ixx, self.Iyy, self.Izz)

        self.rigid_body = RigidBody('rigid_body1', self.P1, self.A,
                                    self.mass, (self.inertia, self.P1))

        self.global_frame1 = VisualizationFrame('global_frame1', self.A,
                                                self.P1, self.shape1)

        self.global_frame2 = VisualizationFrame('global_frame2', self.B,
                                                self.P2, self.shape2)

        self.scene1 = Scene(self.inertial_ref_frame, self.origin,
                            (self.global_frame1, self.global_frame2),
                            name='scene')

        self.particle = Particle('particle1', self.P1, self.mass)

        q = self.coords
        c = cos
        s = sin
        # Here is the dragon ..
        self.transformation_matrix = [
            [c(q[1])*c(q[2]), s(q[2])*c(q[1]), -s(q[1]), 0],
            [s(q[0])*s(q[1])*c(q[2]) - s(q[2])*c(q[0]),
             s(q[0])*s(q[1])*s(q[2]) + c(q[0])*c(q[2]),
             s(q[0])*c(q[1]), 0],
            [s(q[0])*s(q[2]) + s(q[1])*c(q[0])*c(q[2]),
             -s(q[0])*c(q[2]) + s(q[1])*s(q[2])*c(q[0]),
             c(q[0])*c(q[1]), 0],
            [10, 10, 10, 1]]

    def test_vframe_with_rframe(self):
        self.frame1 = VisualizationFrame('frame1', self.inertial_ref_frame,
                                         self.origin, self.shape1)

        assert self.frame1.name == 'frame1'
        assert self.frame1.reference_frame == self.inertial_ref_frame
        assert self.frame1.origin == self.origin
        assert self.frame1.shape is self.shape1

        self.frame1.name = 'frame1_'
        assert self.frame1.name == 'frame1_'

        self.frame1.reference_frame = self.A
        assert self.frame1.reference_frame == self.A

        self.frame1.origin = self.P1
        assert self.frame1.origin == self.P1

        self.frame1.shape = self.shape2
        assert self.frame1.shape is self.shape2

        tr_mat = self.frame1.generate_transformation_matrix(
            self.inertial_ref_frame, self.origin).tolist()
        assert tr_mat == self.transformation_matrix

    def test_vframe_with_rbody(self):
        self.frame2 = VisualizationFrame('frame2', self.rigid_body,
                                         self.shape1)

        assert self.frame2.name == 'frame2'
        assert self.frame2.reference_frame == self.A
        assert self.frame2.origin == self.P1
        assert self.frame2.shape == self.shape1

        self.frame2.name = 'frame2_'
        assert self.frame2.name == 'frame2_'

        self.frame2.reference_frame = self.B
        assert self.frame2.reference_frame == self.B

        self.frame2.origin = self.P2
        assert self.frame2.origin == self.P2

        self.frame2.shape = self.shape2
        assert self.frame2.shape is self.shape2

        self.frame2.reference_frame = self.A
        self.frame2.origin = self.P1
        tr_mat = self.frame2.generate_transformation_matrix(
            self.inertial_ref_frame, self.origin).tolist()
        assert tr_mat == self.transformation_matrix

    def test_vframe_with_particle(self):

        self.frame3 = VisualizationFrame('frame3', self.A, self.particle,
                                         self.shape1)

        assert self.frame3.name == 'frame3'
        assert self.frame3.reference_frame == self.A
        assert self.frame3.origin == self.P1
        assert self.frame3.shape is self.shape1

        self.frame3.name = 'frame3_'
        assert self.frame3.name == 'frame3_'

        self.frame3.reference_frame = self.B
        assert self.frame3.reference_frame == self.B

        self.frame3.origin = self.P2
        assert self.frame3.origin == self.P2

        self.frame3.shape = self.shape2
        assert self.frame3.shape is self.shape2

        self.frame3.reference_frame = self.A
        self.frame3.origin = self.P1
        tr_mat = self.frame3.generate_transformation_matrix(
            self.inertial_ref_frame, self.origin).tolist()
        assert tr_mat == self.transformation_matrix

    def test_vframe_without_name(self):

        self.frame4 = VisualizationFrame(self.inertial_ref_frame, self.origin,
                                         self.shape1)

        assert self.frame4.name == 'unnamed'
        # To check if referenceframe and origin are defined properly without
        # name arg
        assert self.frame4.reference_frame == self.inertial_ref_frame
        assert self.frame4.origin == self.origin
        assert self.frame4.shape is self.shape1

        self.frame4.name = 'frame1_'
        assert self.frame4.name == 'frame1_'

    def test_numeric_transform(self):

        self.list1 = [[0.5000000000000001,
                       0.5,
                       -0.7071067811865475,
                       0.0,
                       -0.14644660940672627,
                       0.8535533905932737,
                       0.5,
                       0.0,
                       0.8535533905932737,
                       -0.14644660940672627,
                       0.5000000000000001,
                       0.0,
                       10.0,
                       10.0,
                       10.0,
                       1.0]]

        self.list2 = [[-0.11518993731879767,
                       0.8178227645734215,
                       -0.563823734943801,
                       0.0,
                       0.1332055011661179,
                       0.5751927992738988,
                       0.8070994598700584,
                       0.0,
                       0.984371663956036,
                       0.017865313009926137,
                       -0.17519491371464685,
                       0.0,
                       20.0,
                       20.0,
                       20.0,
                       1.0]]

        times = [1.0, 2.0, 3.0, 4.0]

        self.global_frame1.generate_transformation_matrix(
            self.inertial_ref_frame, self.origin)
        self.global_frame1.generate_numeric_transform_function(
            self.dynamic, self.parameters)

        assert_allclose(self.global_frame1.evaluate_transformation_matrix(
            self.states, self.param_vals), self.list1)

        self.global_frame2.generate_transformation_matrix(
            self.inertial_ref_frame, self.origin)
        self.global_frame2.generate_numeric_transform_function(
            self.dynamic, self.parameters)

        assert_allclose(self.global_frame2.evaluate_transformation_matrix(
            self.states, self.param_vals), self.list2)

        try:
            import pythreejs
        except ImportError:
            print('pythreejs not installed, skipping related tests.')
        else:
            del pythreejs
            self.global_frame2._create_keyframetrack(
                times, self.states, self.param_vals,
                constant_map=dict(zip(self.parameters, self.param_vals)))
            self.global_frame2._track.name == 'scene/booger.matrix'
            assert_allclose(self.global_frame2._track.times, times)

    def test_perspective_camera(self):

        # Camera is a subclass of VisualizationFrame, but without any specific
        # shape attached. We supply only ReferenceFrame,Point to camera. and it
        # inherits methods from VisualizationFrame

        camera = PerspectiveCamera('camera', self.rigid_body, fov=45, near=1,
                                   far=1000)

        assert camera.name == 'camera'
        assert camera.reference_frame == self.A
        assert camera.origin == self.P1
        assert camera.fov == 45
        assert camera.near == 1
        assert camera.far == 1000

        camera = PerspectiveCamera('camera', self.inertial_ref_frame,
                                   self.particle, fov=45, near=1, far=1000)

        assert camera.name == 'camera'
        assert camera.reference_frame == self.inertial_ref_frame
        assert camera.origin == self.P1
        assert camera.fov == 45
        assert camera.near == 1
        assert camera.far == 1000

        camera = PerspectiveCamera('camera', self.inertial_ref_frame,
                                   self.origin, fov=45, near=1, far=1000)

        assert camera.name == 'camera'
        assert camera.reference_frame == self.inertial_ref_frame
        assert camera.origin == self.origin
        assert camera.fov == 45
        assert camera.near == 1
        assert camera.far == 1000

        camera.name = 'camera1'
        assert camera.name == 'camera1'
        assert camera.__str__() == 'PerspectiveCamera: camera1'
        assert camera.__repr__() == 'PerspectiveCamera'

        camera.reference_frame = self.A
        assert camera.reference_frame == self.A

        camera.origin = self.P1
        assert camera.origin == self.P1

        camera.fov = 30
        assert camera.fov == 30

        camera.near = 10
        assert camera.near == 10

        camera.far = 500
        assert camera.far == 500

        camera = PerspectiveCamera('camera', self.camera_ref_frame,
                                   self.origin, fov=45, near=1, far=1000)
        camera.generate_transformation_matrix(self.inertial_ref_frame,
                                              self.origin)
        camera.generate_numeric_transform_function(self.dynamic,
                                                   self.parameters)
        viz_matrix = camera.evaluate_transformation_matrix(self.states,
                                                           self.param_vals)
        expected_viz_matrix = [[1., 0., 0., 0.,
                                0., 0., 1., 0.,
                                0., -1., 0., 0.,
                                0., 0., 0., 1.]]
        assert_allclose(viz_matrix, expected_viz_matrix)

        scene_dict = camera.generate_scene_dict()

        assert scene_dict[id(camera)]['type'] == 'PerspectiveCamera'
        assert scene_dict[id(camera)]['fov'] == 45
        assert scene_dict[id(camera)]['near'] == 1
        assert scene_dict[id(camera)]['far'] == 1000
        assert scene_dict[id(camera)]['name'] == 'camera'
        assert scene_dict[id(camera)]['init_orientation'] == \
            [1.,  0.,  0.,  0.,
             0.,  0.,  1.,  0.,
             0., -1.,  0.,  0.,
             0.,  0.,  0.,  1.]

        camera1 = PerspectiveCamera(self.inertial_ref_frame, self.origin)
        assert camera1.name == 'unnamed'
        assert camera1.reference_frame == self.inertial_ref_frame
        assert camera1.origin == self.origin
        assert camera1.fov == 45
        assert camera1.near == 1
        assert camera1.far == 1000

    def test_orthographic_camera(self):

        # As compared to Perspective Camera, Orthographic Camera doesn't have
        # fov, instead the left,right,top and bottom faces are adjusted by the
        # Scene width, and height

        # Testing with rigid_body
        camera = OrthoGraphicCamera('camera', self.rigid_body, near=1,
                                    far=1000)

        assert camera.name == 'camera'
        assert camera.reference_frame == self.A
        assert camera.origin == self.P1
        assert camera.near == 1
        assert camera.far == 1000

        # Testing with reference_frame, particle
        camera = OrthoGraphicCamera('camera', self.inertial_ref_frame,
                                    self.particle, near=1, far=1000)

        assert camera.name == 'camera'
        assert camera.reference_frame == self.inertial_ref_frame
        assert camera.origin == self.P1
        assert camera.near == 1
        assert camera.far == 1000

        # Testing with reference_frame, point ...
        camera = OrthoGraphicCamera('camera', self.inertial_ref_frame,
                                    self.origin, near=1, far=1000)

        assert camera.name == 'camera'
        assert camera.reference_frame == self.inertial_ref_frame
        assert camera.origin == self.origin
        assert camera.near == 1
        assert camera.far == 1000

        camera.name = 'camera1'
        assert camera.name == 'camera1'
        assert camera.__str__() == 'OrthoGraphicCamera: camera1'
        assert camera.__repr__() == 'OrthoGraphicCamera'

        camera.reference_frame = self.A
        assert camera.reference_frame == self.A

        camera.origin = self.P1
        assert camera.origin == self.P1

        camera.near = 10
        assert camera.near == 10

        camera.far = 500
        assert camera.far == 500

        camera1 = OrthoGraphicCamera(self.inertial_ref_frame, self.origin)
        assert camera1.name == 'unnamed'
        assert camera1.reference_frame == self.inertial_ref_frame
        assert camera1.origin == self.origin
        assert camera1.near == 1
        assert camera1.far == 1000

        camera = OrthoGraphicCamera('camera', self.camera_ref_frame,
                                    self.origin, near=1, far=1000)
        camera.generate_transformation_matrix(self.inertial_ref_frame,
                                              self.origin)
        camera.generate_numeric_transform_function(self.dynamic,
                                                   self.parameters)
        viz_matrix = camera.evaluate_transformation_matrix(self.states,
                                                           self.param_vals)
        expected_viz_matrix = [[1., 0., 0., 0.,
                                0., 0., 1., 0.,
                                0., -1., 0., 0.,
                                0., 0., 0., 1.]]
        assert_allclose(viz_matrix, expected_viz_matrix)

        scene_dict = camera.generate_scene_dict()

        assert scene_dict[id(camera)]['type'] == 'OrthoGraphicCamera'
        assert scene_dict[id(camera)]['near'] == 1
        assert scene_dict[id(camera)]['far'] == 1000
        assert scene_dict[id(camera)]['name'] == 'camera'
        assert scene_dict[id(camera)]['init_orientation'] == \
            [1.,  0.,  0.,  0.,
             0.,  0.,  1.,  0.,
             0., -1.,  0.,  0.,
             0.,  0.,  0.,  1.]

    def test_point_light(self):

        # Testing with rigid-body ..
        light = PointLight('light', self.rigid_body, color='blue')

        assert light.name == 'light'
        assert light.reference_frame == self.A
        assert light.origin == self.P1
        assert light.color == 'blue'

        # Testing with reference_frame, particle ..
        light = PointLight('light', self.inertial_ref_frame, self.particle,
                           color='blue')

        assert light.name == 'light'
        assert light.reference_frame == self.inertial_ref_frame
        assert light.origin == self.P1
        assert light.color == 'blue'

        # Testing with reference_frame, point ..
        light = PointLight('light', self.inertial_ref_frame, self.origin,
                           color='blue')

        assert light.name == 'light'
        assert light.reference_frame == self.inertial_ref_frame
        assert light.origin == self.origin
        assert light.color == 'blue'

        light.name = 'light1'
        assert light.name == 'light1'
        assert light.__str__() == 'PointLight: light1'
        assert light.__repr__() == 'PointLight'

        light.reference_frame = self.A
        assert light.reference_frame == self.A

        light.origin = self.P1
        assert light.origin == self.P1

        light.color = 'red'
        assert light.color == 'red'

        # Test unnamed
        light1 = PointLight(self.inertial_ref_frame, self.origin)
        assert light1.name == 'unnamed'
        assert light1.reference_frame == self.inertial_ref_frame
        assert light1.origin == self.origin
        assert light1.color == 'white'

    def test_scene_init(self):

        self.scene2 = Scene(self.inertial_ref_frame, self.origin,
                            self.global_frame1, self.global_frame2,
                            name='scene')

        assert self.scene2.name == 'scene'
        assert self.scene2.reference_frame == self.inertial_ref_frame
        assert self.scene2.origin == self.origin
        assert self.scene2.visualization_frames[0] is self.global_frame1
        assert self.scene2.visualization_frames[1] is self.global_frame2

        self.scene2.name = 'scene1'
        assert self.scene2.name == 'scene1'

        self.scene2.reference_frame = self.A
        assert self.scene2.reference_frame == self.A
