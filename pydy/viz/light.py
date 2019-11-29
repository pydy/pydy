from sympy.matrices.expressions import Identity

from .visualization_frame import VisualizationFrame
from ..utils import sympy_equal_to_or_newer_than

__all__ = ['PointLight']

class PointLight(VisualizationFrame):
    """
    Creates a PointLight for the visualization
    The PointLight is inherited from VisualizationFrame,

    It can also be attached to dynamics objects, hence we can
    get a moving Light. All the transformation matrix generation
    methods are applicable to a PointLight.
    Like VisualizationFrame,
    It can also be initialized using:
    1)Rigidbody
    2)ReferenceFrame, Point
    3)ReferenceFrame, Particle
    Either one of these must be supplied during initialization

    Unlike VisualizationFrame, It doesnt require a Shape argument.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialises a PointLight object.
        To initialize a point light, we need to supply
        a name(optional), a reference frame, and a point.


        Parameters
        ==========
        name : str, optional
            Name assigned to VisualizationFrame, default is unnamed
        reference_frame : ReferenceFrame
            A reference_frame with respect to which all orientations of the
            shape takes place, during visualizations/animations.
        origin : Point
            A point with respect to which all the translations of the shape
            takes place, during visualizations/animations.
        rigidbody : RigidBody
            A rigidbody whose reference frame and mass center are to be
            assigned as reference_frame and origin of the
            VisualizationFrame.
        particle : Particle
            A particle whose point is assigned as origin of the
            VisualizationFrame.

        Examples
        ========
        >>> from pydy.viz import PointLight
        >>> from sympy.physics.mechanics import \
                               ReferenceFrame, Point, RigidBody, \
                                Particle, inertia
        >>> from sympy import symbols
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')
        >>> #initializing with reference frame, point
        >>> light = PointLight('light', I, O)
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> # Initializing with a rigidbody ..
        >>> light = PointLight('frame2', rbody)
        >>> Pa = Particle('Pa', O, mass)
        >>> #initializing with Particle, reference_frame ...
        >>> light = PointLight('frame3', I, Pa)
        """

        try:
            self._color = kwargs['color']
        except KeyError:
            self._color = 'white'

        #Now we use same approach as in VisualizationFrame
        #for setting reference_frame and origin
        i = 0
        #If first arg is not str, name the visualization frame 'unnamed'
        if isinstance(args[i], str):
            self._name = args[i]
            i += 1
        else:
            self._name = 'unnamed'

        try:
            if sympy_equal_to_or_newer_than('1.0'):
                self._reference_frame = args[i].frame
            else:
                self._reference_frame = args[i].get_frame()

            self._origin = args[i].masscenter

        except AttributeError:
            #It is not a rigidbody, hence this arg should be a
            #reference frame
            try:
                dcm = args[i]._dcm_dict
                self._reference_frame = args[i]
                i += 1
            except AttributeError:
                raise TypeError(''' A ReferenceFrame is to be supplied
                                   before a Particle/Point. ''')

            #Now next arg can either be a Particle or point
            try:
                self._origin = args[i].point
            except AttributeError:
                self._origin = args[i]

        #basic thing required, transform matrix
        self._transform = Identity(4).as_mutable()

    def __str__(self):
        return 'PointLight: ' + self._name

    def __repr__(self):
        return 'PointLight'

    @property
    def color(self):
        """
        Color of Light.
        """
        return self._color

    @color.setter
    def color(self, new_color):
        """
        Sets the color of Light.

        """
        if not isinstance(new_color, str):
            raise TypeError('''Color should be a valid str object''')
        else:
            self._color = new_color

    def color_in_rgb(self):
        """
        Returns the rgb value of the
        defined light color.
        """
        return self._color_rgb

    def generate_scene_dict(self):
        """
        This method generates information for a static
        visualization in the initial conditions, in the form
        of dictionary. This contains light parameters followed
        by an init_orientation Key.

        Before calling this method, all the transformation matrix
        generation methods should be called, or it will give an error.
        Returns
        =======
        A dict with following Keys:

        1. name: name for the camera
        2. color: Color of the light
        3. init_orientation: Initial orientation
           of the light object

        """
        scene_dict = { id(self): {} }
        scene_dict[id(self)]['name'] = self.name
        scene_dict[id(self)]['type'] = self.__repr__()
        scene_dict[id(self)]['color'] = self.color
        scene_dict[id(self)]["simulation_id"] = id(self)
        scene_dict[id(self)]["init_orientation"] = self._visualization_matrix[0]

        return scene_dict

    def generate_simulation_dict(self):
        """
        Generates the simulation information for this Light object.
        It maps the simulation data information to the
        scene information via a unique id.

        Before calling this method, all the transformation matrix
        generation methods should be called, or it will give an error.

        Returns
        =======

        A dictionary containing list of 4x4 matrices mapped to
        the unique id as the key.

        """
        simulation_dict = {}
        try:
            simulation_dict[id(self)] = self._visualization_matrix

        except:
            raise RuntimeError("Please call the numerical ",
                               "transformation methods, ",
                               "before generating visualization dict")


        return simulation_dict
