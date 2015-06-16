import warnings
from sympy.matrices.expressions import Identity

from ..utils import PyDyUserWarning
from .visualization_frame import VisualizationFrame

__all__ = ['PerspectiveCamera', 'OrthoGraphicCamera']

warnings.simplefilter('once', PyDyUserWarning)


class PerspectiveCamera(VisualizationFrame):
    """Creates a Perspective Camera for visualization. The camera is
    inherited from VisualizationFrame. It can be attached to dynamics
    objects, hence we can get a moving camera. All the transformation matrix
    generation methods are applicable to a Perspective Camera.

    Like VisualizationFrame, it can also be initialized using:

       1. Rigidbody
       2. ReferenceFrame, Point
       3. ReferenceFrame, Particle

    Either one of these must be supplied during initialization. Unlike
    VisualizationFrame, it does not require a Shape argument.

    Parameters
    ==========
    name : str
        A name for the PerspectiveCamera(optional). Default is 'unnamed'
    fov : float, default=45.0
        Field Of View, It determines the angle between the top and bottom of
        the viewable area (in degrees).
    near : float
        The distance of near plane of the PerspectiveCamera. All objects
        closer to this distance are not displayed.
    far : int or float
        The distance of far plane of the PerspectiveCamera. All objects
        farther than this distance are not displayed.

    """

    def __init__(self, *args, **kwargs):
        """Initialises a PerspectiveCamera object. To initialize a
        visualization frame, one needs to supply a name (optional), a
        reference frame, a point, field of view (fov) (optional), near plane
        distance (optional) and far plane distance (optional).

        Examples
        ========
        >>> from pydy.viz import VisualizationFrame, Shape
        >>> from sympy.physics.mechanics import (ReferenceFrame, Point,
        ...                                      RigidBody, Particle,
        ...                                      inertia)
        >>> from sympy import symbols
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')
        >>> shape = Shape()
        >>> # initializing with reference frame, point
        >>> camera1 = PerspectiveCamera('frame1', I, O)
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> # Initializing with a rigidbody ..
        >>> camera2 = PerspectiveCamera('frame2', rbody)
        >>> Pa = Particle('Pa', O, mass)
        >>> # initializing with Particle, reference_frame ...
        >>> camera3 = PerspectiveCamera('frame3', I, Pa)
        """

        msg = ("Rotation of Perspective Camera does not work "
               "properly in the visualiser.")
        warnings.warn(msg, PyDyUserWarning)

        try:
            self._fov = kwargs['fov']
        except KeyError:
            self._fov = 45.0

        try:
            self._near = kwargs['near']
        except KeyError:
            self._near = 1.0

        try:
            self._far = kwargs['far']
        except KeyError:
            self._far = 1000.0

        # Now we use same approach as in VisualizationFrame for setting
        # reference_frame and origin
        i = 0
        # If first arg is not str, name the visualization frame 'unnamed'
        if isinstance(args[i], str):
            self._name = args[i]
            i += 1
        else:
            self._name = 'unnamed'

        try:
            self._reference_frame = args[i].get_frame()
            self._origin = args[i].get_masscenter()

        except AttributeError:
            # It is not a rigidbody, hence this arg should be a reference
            # frame
            try:
                # TODO : dcm is never used.
                dcm = args[i]._dcm_dict
                self._reference_frame = args[i]
                i += 1
            except AttributeError:
                raise TypeError('A ReferenceFrame is to be supplied '
                                'before a Particle/Point.')

            # Now next arg can either be a Particle or point
            try:
                self._origin = args[i].get_point()
            except AttributeError:
                self._origin = args[i]

        # basic thing required, transform matrix
        self._transform = Identity(4).as_mutable()

    def __str__(self):
        return 'PerspectiveCamera: ' + self._name

    def __repr__(self):
        return 'PerspectiveCamera'

    @property
    def fov(self):
        """
        attribute for Field Of view of a PerspectiveCamera.
        Default value is 45 degrees
        """
        return self._fov

    @fov.setter
    def fov(self, new_fov):
        if not isinstance(new_fov, (int, str)):
            raise TypeError('fov should be supplied in int or float')
        else:
            self._fov = new_fov

    @property
    def near(self):
        """
        attribute for Near Plane distance of a PerspectiveCamera.
        Default value is 1
        """
        return self._near

    @near.setter
    def near(self, new_near):
        if not isinstance(new_near, (int, str)):
            raise TypeError('near should be supplied in int or float')
        else:
            self._near = new_near

    @property
    def far(self):
        """
        Attribute for Far Plane distance of a PerspectiveCamera. The default
        value is ``1000.0``.
        """
        return self._far

    @far.setter
    def far(self, new_far):
        if not isinstance(new_far, (int, str)):
            raise TypeError('far should be supplied in int or float')
        else:
            self._far = new_far

    def generate_scene_dict(self):
        """This method generates information for a static visualization in
        the initial conditions, in the form of dictionary. This contains
        camera parameters followed by an init_orientation Key.

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
        scene_dict = {id(self): {}}
        scene_dict[id(self)]['name'] = self.name
        scene_dict[id(self)]['type'] = self.__repr__()
        scene_dict[id(self)]['fov'] = self.fov
        scene_dict[id(self)]['near'] = self.near
        scene_dict[id(self)]['far'] = self.far
        scene_dict[id(self)]["simulation_id"] = id(self)
        scene_dict[id(self)]["init_orientation"] = self._visualization_matrix[0]

        return scene_dict

    def generate_simulation_dict(self):
        """Generates the simulation information for this visualization
        frame.  It maps the simulation data information to the scene
        information via a unique id.

        Before calling this method, all the transformation matrix generation
        methods should be called, or it will give an error.

        Returns
        =======
        simulation_dict : dictionary
            A dictionary containing list of 4x4 matrices mapped to the
            unique id as the key.

        """
        simulation_dict = {}
        try:
            simulation_dict[id(self)] = self._visualization_matrix

        except:
            raise RuntimeError("Please call the numerical "
                               "transformation methods, "
                               "before generating visualization dict.")

        return simulation_dict


class OrthoGraphicCamera(VisualizationFrame):
    """Creates a OrthoGraphic Camera for visualization. The camera is
    inherited from ``VisualizationFrame``. It can be attached to dynamics
    objects, hence we can get a moving camera. All the transformation matrix
    generation methods are applicable to a Perspective Camera.

    Like VisualizationFrame, it can also be initialized using:

       1. :role:`Rigidbody`
       2. ReferenceFrame, Point
       3. ReferenceFrame, Particle

    Either one of these must be supplied during initialization. Unlike
    VisualizationFrame, it doesnt require a Shape argument.

    Parameters
    ==========
    name : str, optional, default='unnamed'
        A name for the PerspectiveCamera.
    near : float, optional, default=
        The distance of near plane of the PerspectiveCamera. All objects
        closer to this distance are not displayed.
    far : float, optional, default=
        The distance of far plane of the PerspectiveCamera. All objects
        farther than this distance are not displayed.

    """

    def __init__(self, *args, **kwargs):
        """
        Initialises an OrthoGraphicCamera object. To initialize a
        visualization frame, we need to supply a name (optional), a
        reference frame, a point, near plane distance (optional) and far
        plane distance (optional).

        Examples
        ========
        >>> from pydy.viz import OrthoGraphicCamera
        >>> from sympy.physics.mechanics import (ReferenceFrame, Point,
        ...                                      RigidBody, Particle,
        ...                                      inertia)
        >>> from sympy import symbols
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')
        >>> shape = Shape()
        >>> # Initializing with ReferenceFrame, Point
        >>> camera1 = OrthoGraphicCamera('frame1', I, O)
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> # Initializing with a Rigidbody
        >>> camera2 = OrthoGraphicCamera('frame2', rbody)
        >>> Pa = Particle('Pa', O, mass)
        >>> # Initializing with Particle, ReferenceFrame
        >>> camera3 = OrthoGraphicCamera('frame3', I, Pa)
        """
        try:
            self._near = kwargs['near']
        except KeyError:
            self._near = 1

        try:
            self._far = kwargs['far']
        except KeyError:
            self._far = 1000

        # Now we use same approach as in VisualizationFrame for setting
        # reference_frame and origin
        i = 0
        # If first arg is not str, name the visualization frame 'unnamed'
        if isinstance(args[i], str):
            self._name = args[i]
            i += 1
        else:
            self._name = 'unnamed'

        try:
            self._reference_frame = args[i].get_frame()
            self._origin = args[i].get_masscenter()

        except AttributeError:
            # It is not a rigidbody, hence this arg should be a reference
            # frame.
            self._reference_frame = args[i]
            i += 1

            # Now next arg can either be a Particle or point
            try:
                self._origin = args[i].get_point()
            except AttributeError:
                self._origin = args[i]

        # basic thing required, transform matrix
        self._transform = Identity(4).as_mutable()

    def __str__(self):
        return 'OrthoGraphicCamera: ' + self._name

    def __repr__(self):
        return 'OrthoGraphicCamera'

    @property
    def near(self):
        """Attribute for Near Plane distance of an OrthoGraphicCamera.
        Default value is 1.0
        """
        return self._near

    @near.setter
    def near(self, new_near):
        if not isinstance(new_near, (int, str)):
            raise TypeError('near should be supplied in int or float')
        else:
            self._near = new_near

    @property
    def far(self):
        """Attribute for Far Plane distance of an OrthoGraphicCamera.
        Default value is 1000.0.
        """
        return self._far

    @far.setter
    def far(self, new_far):
        if not isinstance(new_far, (int, str)):
            raise TypeError('far should be supplied in int or float')
        else:
            self._far = new_far

    def generate_scene_dict(self):
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
        scene_dict = {id(self): {}}
        scene_dict[id(self)]['name'] = self.name
        scene_dict[id(self)]['type'] = self.__repr__()
        scene_dict[id(self)]['near'] = self.near
        scene_dict[id(self)]['far'] = self.far
        scene_dict[id(self)]["simulation_id"] = id(self)
        scene_dict[id(self)]["init_orientation"] = self._visualization_matrix[0]

        return scene_dict

    def generate_simulation_dict(self):
        """Generates the simulation information for this visualization
        frame.  It maps the simulation data information to the scene
        information via a unique id.

        Returns
        =======

        A dictionary containing list of 4x4 matrices mapped to the unique id
        as the key.

        """
        simulation_dict = {}
        try:
            simulation_dict[id(self)] = self._visualization_matrix

        except:
            raise RuntimeError("Please call the numerical ",
                               "transformation methods, ",
                               "before generating visualization dict")

        return simulation_dict
