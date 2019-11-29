#!/usr/bin/env python

# standard lib
import warnings

# local
from ..utils import PyDyUserWarning
from .shapes import Shape
from .visualization_frame import VisualizationFrame

__all__ = ['PerspectiveCamera', 'OrthoGraphicCamera']

warnings.simplefilter('once', PyDyUserWarning)


class PerspectiveCamera(VisualizationFrame):
    """Creates a perspective camera for use in a scene. The camera is inherited
    from ``VisualizationFrame``, and thus behaves similarly. It can be attached
    to dynamics objects, hence we can get a moving camera. All the
    transformation matrix generation methods are applicable to a
    ``PerspectiveCamera``.

    """

    def __init__(self, *args, **kwargs):
        """Initialises a PerspectiveCamera object. To initialize a
        PerspectiveCamera, one needs to supply a name (optional), a reference
        frame, a point, field of view (fov) (optional), near plane distance
        (optional) and far plane distance (optional).

        Like ``VisualizationFrame``, it can also be initialized using one of
        these three argument sequences:

        Rigidbody
        ``PerspectiveCamera(rigid_body)``
        ReferenceFrame, Point
        ``PerspectiveCamera(ref_frame, point)``
        ReferenceFrame, Particle
        ``PerspectiveCamera(ref_frame, particle)``

        Note that you can also supply and optional name as the first positional
        argument, e.g.::

           ``PerspectiveCamera('camera_name', rigid_body)``

        Additional optional keyword arguments are below:

        Parameters
        ==========
        fov : float, default=45.0
            Field Of View, It determines the angle between the top and bottom
            of the viewable area (in degrees).
        near : float
            The distance of near plane of the PerspectiveCamera. All objects
            closer to this distance are not displayed.
        far : int or float
            The distance of far plane of the PerspectiveCamera. All objects
            farther than this distance are not displayed.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import (ReferenceFrame, Point,
        ...                                      RigidBody, Particle,
        ...                                      inertia)
        >>> from pydy.viz import PerspectiveCamera
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')

        >>> # initializing with reference frame, point
        >>> camera1 = PerspectiveCamera('frame1', I, O)

        >>> # Initializing with a RigidBody
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> camera2 = PerspectiveCamera('frame2', rbody)

        >>> # initializing with Particle, reference_frame
        >>> Pa = Particle('Pa', O, mass)
        >>> camera3 = PerspectiveCamera('frame3', I, Pa)

        """

        # NOTE: This allows us to use inhertiance even though cameras don't
        # need a shape. In the future, this could be a camera shape that could
        # be made visible in the scene (only important for multiple cameras).
        args = list(args) + [Shape()]

        super(PerspectiveCamera, self).__init__(*args)

        self.fov = 45.0
        self.near = 1.0
        self.far = 1000.0

        for k, v in kwargs.items():
            setattr(self, k, v)

    def __str__(self):
        return 'PerspectiveCamera: ' + self.name

    def __repr__(self):
        return 'PerspectiveCamera'

    @property
    def fov(self):
        return self._fov

    @fov.setter
    def fov(self, new_fov):
        self._fov = float(new_fov)

    @property
    def near(self):
        return self._near

    @near.setter
    def near(self, new_near):
        self._near = float(new_near)

    @property
    def far(self):
        return self._far

    @far.setter
    def far(self, new_far):
        self._far = float(new_far)

    def generate_scene_dict(self, **kwargs):
        """This method generates information for a static visualization in the
        initial conditions, in the form of dictionary. This contains camera
        parameters followed by an init_orientation key.

        Before calling this method, all the transformation matrix generation
        methods should be called, or it will give an error.

        Returns
        =======
        A dict with following Keys:

        1. name: name for the camera
        2. fov: Field of View value of the camera
        3. near: near value of the camera
        4. far: far value of the camera
        5. init_orientation: Initial orientation of the camera

        """
        scene_dict = super(PerspectiveCamera, self).generate_scene_dict(**kwargs)
        scene_dict[id(self)]['type'] = self.__class__.__name__
        scene_dict[id(self)]['fov'] = self.fov
        scene_dict[id(self)]['near'] = self.near
        scene_dict[id(self)]['far'] = self.far

        return scene_dict


class OrthoGraphicCamera(VisualizationFrame):
    """Creates a orthographic camera for use in a scene. The camera is
    inherited from ``VisualizationFrame``, and thus behaves similarly. It can
    be attached to dynamics objects, hence we can get a moving camera. All the
    transformation matrix generation methods are applicable to a
    ``OrthoGraphicCameraCamera``.


    """

    def __init__(self, *args, **kwargs):
        """Initialises a OrthoGraphicCameraCamera object. To initialize a
        OrthoGraphicCameraCamera, one needs to supply a name (optional), a
        reference frame, a point, field of view (fov) (optional), near plane
        distance (optional) and far plane distance (optional).

        Like ``VisualizationFrame``, it can also be initialized using one of
        these three argument sequences:

        Rigidbody
           ``OrthoGraphicCameraCamera(rigid_body)``
        ReferenceFrame, Point
           ``OrthoGraphicCameraCamera(ref_frame, point)``
        ReferenceFrame, Particle
           ``OrthoGraphicCameraCamera(ref_frame, particle)``

        Note that you can also supply and optional name as the first positional
        argument, e.g.::

            OrthoGraphicCameraCamera('camera_name', rigid_body)

        Additional optional keyword arguments are below:

        Parameters
        ==========
        near : float
            The distance of near plane of the OrthoGraphicCameraCamera. All
            objects closer to this distance are not displayed.
        far : int or float
            The distance of far plane of the OrthoGraphicCameraCamera. All
            objects farther than this distance are not displayed.

        Examples
        ========

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import (ReferenceFrame, Point,
        ...                                      RigidBody, Particle,
        ...                                      inertia)
        >>> from pydy.viz import OrthoGraphicCameraCamera
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')

        >>> # initializing with reference frame, point
        >>> camera1 = OrthoGraphicCameraCamera('frame1', I, O)

        >>> # Initializing with a RigidBody
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> camera2 = OrthoGraphicCameraCamera('frame2', rbody)

        >>> # initializing with Particle, reference_frame
        >>> Pa = Particle('Pa', O, mass)
        >>> camera3 = OrthoGraphicCameraCamera('frame3', I, Pa)

        """
        # NOTE: This allows us to use inhertiance even though cameras don't
        # need a shape. In the future, this could be a camera shape that could
        # be made visible in the scene (only important for multiple cameras).
        args = list(args) + [Shape()]

        super(OrthoGraphicCamera, self).__init__(*args)

        self.near = 1.0
        self.far = 1000.0

        for k, v in kwargs.items():
            setattr(self, k, v)

    def __str__(self):
        return self.__class__.__name__ + ': ' + self.name

    def __repr__(self):
        return self.__class__.__name__

    @property
    def near(self):
        return self._near

    @near.setter
    def near(self, new_near):
        self._near = float(new_near)

    @property
    def far(self):
        return self._far

    @far.setter
    def far(self, new_far):
        self._far = float(new_far)

    def generate_scene_dict(self, **kwargs):
        """
        This method generates information for a static visualization in the
        initial conditions, in the form of dictionary. This contains camera
        parameters followed by an init_orientation Key.

        Returns
        =======
        scene_dict : dictionary
           A dict with following Keys:

              1. name: name for the camera
              2. near: near value of the camera
              3. far: far value of the camera
              4. init_orientation: Initial orientation of the camera

        """
        scene_dict = super(OrthoGraphicCamera, self).generate_scene_dict(**kwargs)
        scene_dict[id(self)]['type'] = self.__class__.__name__
        scene_dict[id(self)]['near'] = self.near
        scene_dict[id(self)]['far'] = self.far

        return scene_dict
