__all__ = ['VisualizationFrame']

import sys
from collections.abc import Iterator
import numpy as np
from sympy import Dummy, lambdify
from sympy.matrices.expressions import Identity
from sympy.physics.mechanics import Point, ReferenceFrame
try:
    import pythreejs as p3js
except ImportError:
    p3js = None

from .shapes import Shape
from ..utils import sympy_equal_to_or_newer_than


class VisualizationFrame(object):
    """
    A VisualizationFrame represents an object that you want to visualize.
    It allows you to easily associate a reference frame and a point
    with a shape.

    A VisualizationFrame can be attached to only one Shape Object.
    It can be nested, i.e we can add/remove multiple visualization frames to
    one visualization frame. On adding the parent frame to the
    Scene object, all the children of the parent visualization frame
    are also added, and hence can be visualized and animated.

    A VisualizationFrame needs to have a ReferenceFrame, and a Point
    for it to form transformation matrices for visualization and
    animations.

    The ReferenceFrame and Point are required to be provided during
    initialization. They can be supplied in the form of any one of these:

    1)reference_frame, point argument.
    2)a RigidBody argument
    3)reference_frame, particle argument.

    In addition to these arguments, A shape argument is also required.
    """
    def __init__(self, *args):
        """
        To initialize a visualization frame a ReferenceFrame,
        Point, and Shape are required. These ReferenceFrame
        and Point can be passed provided in three ways:

        1) RigidBody: the RigidBody's frame and mass center are used.
        2) ReferenceFrame and a Particle: The Particle's Point is used.
        3) ReferenceFrame and a Point

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
        shape : Shape
            A shape to be attached to the VisualizationFrame

        Examples
        ========
        >>> from pydy.viz import VisualizationFrame, Sphere
        >>> from sympy.physics.mechanics import \
                               ReferenceFrame, Point, RigidBody, \
                                Particle, inertia
        >>> from sympy import symbols
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')
        >>> shape = Sphere(5)
        >>> #initializing with reference frame, point
        >>> frame1 = VisualizationFrame('frame1', I, O, shape)
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> # Initializing with a rigidbody ..
        >>> frame2 = VisualizationFrame('frame2', rbody, shape)
        >>> Pa = Particle('Pa', O, mass)
        >>> #initializing with Particle, reference_frame ...
        >>> frame3 = VisualizationFrame('frame3', I, Pa, shape)
        """
        #Last arg should be a Shape ..
        if isinstance(args[-1], Shape):
            self._shape = args[-1]
        else:
            raise TypeError("Please provide a valid shape object as the last "
                            " positional argument.")
        i = 0
        # If first arg is not str, name the visualization frame 'unnamed'
        if isinstance(args[i], str):
            self.name = args[i]
            i += 1
        else:
            self.name = 'unnamed'

        try:
            if sympy_equal_to_or_newer_than('1.0'):
                self.reference_frame = args[i].frame
            else:
                self.reference_frame = args[i].get_frame()
            self.origin = args[i].masscenter

        except AttributeError:
            #It is not a rigidbody, hence this arg should be a
            #reference frame
            try:
                dcm = args[i]._dcm_dict
                self.reference_frame = args[i]
                i += 1
            except AttributeError:
                raise TypeError(''' A ReferenceFrame is to be supplied
                                   before a Particle/Point. ''')

            #Now next arg can either be a Particle or point
            try:
                self.origin = args[i].point

            except AttributeError:
                self.origin = args[i]

    #setting attributes ..
    def __str__(self):
        return 'VisualizationFrame ' + self.name

    def __repr__(self):
        return 'VisualizationFrame'

    @property
    def name(self):
        """
        Name of the VisualizationFrame.
        """
        return self._name

    @name.setter
    def name(self, new_name):
        """
        Sets the name of the VisualizationFrame.

        """
        if not isinstance(new_name, str):
            raise TypeError('''Name should be a str object''')
        else:
            self._name = new_name

    @property
    def origin(self):
        """
        Origin of the VisualizationFrame,
        with respect to which all translational transformations
        take place.
        """
        return self._origin

    @origin.setter
    def origin(self, new_origin):
        """
        Sets the origin of the VisualizationFrame.
        """
        if not isinstance(new_origin, Point):
            raise TypeError('''origin should be a valid Point Object''')
        else:
            self._origin = new_origin

    @property
    def reference_frame(self):
        """
        reference_frame of the VisualizationFrame,
        with respect to which all rotational/orientational
        transformations take place.
        """
        return self._reference_frame

    @reference_frame.setter
    def reference_frame(self, new_reference_frame):
        if not isinstance(new_reference_frame, ReferenceFrame):
            raise TypeError('''reference_frame should be a valid
                                ReferenceFrame object.''')
        else:
            self._reference_frame = new_reference_frame

    @property
    def shape(self):
        """
        shape in the VisualizationFrame.
        A shape attached to the visualization frame.
        NOTE: Only one shape can be attached to a visualization frame.
        """
        return self._shape

    @shape.setter
    def shape(self, new_shape):
        """
        Sets the shape for VisualizationFrame.
        """
        if not isinstance(new_shape, Shape):
            raise TypeError('''shape should be a valid Shape object.''')
        else:
            self._shape = new_shape

    def generate_transformation_matrix(self, reference_frame, point):
        """Generates a symbolic transformation matrix, with respect to the
        provided reference frame and point.

        Parameters
        ==========
        reference_frame : ReferenceFrame
            A reference_frame with respect to which transformation matrix is
            generated.
        point : Point
            A point with respect to which transformation matrix is
            generated.

        Returns
        =======
        A 4 x 4 SymPy matrix, containing symbolic expressions describing the
        transformation as a function of time.

        """
        rotation_matrix = self.reference_frame.dcm(reference_frame)
        self._transform = Identity(4).as_mutable()
        self._transform[:3, :3] = rotation_matrix

        point_vector = self.origin.pos_from(point)
        try:
            self._transform[3, :3] = point_vector.to_matrix(reference_frame).T
        except AttributeError:
            # In earlier versions of sympy, 'Vector' object has no attribute
            # 'to_matrix'.
            self._transform[3, 0] = point_vector.dot(reference_frame.x)
            self._transform[3, 1] = point_vector.dot(reference_frame.y)
            self._transform[3, 2] = point_vector.dot(reference_frame.z)
        return self._transform

    def generate_numeric_transform_function(self, dynamic_variables,
                                            constant_variables):
        """Returns a function which can compute the numerical values of the
        transformation matrix given the numerical dynamic variables (i.e.
        functions of time or states) and the numerical system constants.


        Parameters
        ==========
        dynamic_variables : list of sympy.Functions(time)
            All of the dynamic symbols used in defining the orientation and
            position of this visualization frame.
        constant_variables : list of sympy.Symbols
            All of the constants used in defining the orientation and
            position of this visualization frame.

        Returns
        =======
        numeric_transform : list of functions
            A list of functions which return the numerical transformation
            for each element in the transformation matrix.

        """

        dummy_symbols = [Dummy() for i in dynamic_variables]
        dummy_dict = dict(zip(dynamic_variables, dummy_symbols))
        transform = self._transform.subs(dummy_dict).reshape(16, 1)
        dummy_symbols.extend(constant_variables)

        # Create a numeric transformation for each element in the transformation
        # matrix. We cannot lambdify the transformation matrix as calling
        # lambdify of a constant expression returns a scalar, even if the
        # lambdify function arguments are sequences:
        # https://github.com/sympy/sympy/issues/5642
        self._numeric_transform = []
        for i in range(16):
            t = transform[i]
            if t.has(Dummy):
                f = lambdify(dummy_symbols, t, modules='numpy')
            else:
                f = lambdify(constant_variables, t, modules='numpy')
            self._numeric_transform.append(f)
        return self._numeric_transform

    def evaluate_transformation_matrix(self, dynamic_values, constant_values):
        """Returns the numerical transformation matrices for each time step.

        Parameters
        ----------
        dynamic_values : array_like, shape(m,) or shape(n, m)
            The m state values for each n time step.
        constant_values : array_like, shape(p,)
            The p constant parameter values of the system.

        Returns
        -------
        transform_matrix : numpy.array, shape(n, 16)
            A 4 x 4 transformation matrix for each time step.

        """
        #If states is instance of numpy array, well and good.
        #else convert it to one:

        states = np.squeeze(np.array(dynamic_values))
        if not isinstance(constant_values, Iterator):
            constant_values = list(constant_values)
        if len(states.shape) > 1:
            n = states.shape[0]
            args = []
            for a in np.split(states, states.shape[1], 1):
                args.append(np.squeeze(a))
            for a in constant_values:
                args.append(np.repeat(a, n))
        else:
            n = 1
            args = np.hstack((states, constant_values))

        new = np.zeros((n, 16))
        for i, t in enumerate(self._numeric_transform):
            if callable(t):
                try:
                    new[:, i] = t(*args)
                except TypeError:
                    # dynamic values are not necessary so pass only constant
                    # values into transform function
                    new[:, i] = np.repeat(t(*constant_values), n)
            else:
                new[:, i] = np.repeat(t, n)
        self._visualization_matrix = new.tolist()
        return self._visualization_matrix

    def generate_scene_dict(self, constant_map={}):
        """
        This method generates information for a static
        visualization in the initial conditions, in the form
        of dictionary. This contains shape information
        from `Shape.generate_dict()` followed by an
        init_orientation Key.

        Before calling this method, all the transformation matrix
        generation methods should be called, or it will give an error.

        Parameters
        ==========
        constant_map : dictionary
            Constant map is required when Shape contains sympy expressions.This
            dictionary maps sympy expressions/symbols to numerical values(floats)

        Returns
        =======
        A dictionary built with a call to `Shape.generate_dict`.
        Additional keys included in the dict are following:

        1. init_orientation: Specifies the initial orientation
           of the `VisualizationFrame`.
        2. reference_frame_name: Name(str) of the reference_frame
           attached to this VisualizationFrame.
        3. simulation_id: an arbitrary integer to map scene description
           with the simulation data.


        """
        scene_dict = { id(self): {} }
        scene_dict[id(self)] = self.shape.generate_dict(constant_map=constant_map)
        scene_dict[id(self)]['name'] = self.name
        scene_dict[id(self)]["reference_frame_name"] = str(self.reference_frame)
        scene_dict[id(self)]["simulation_id"] = id(self)

        try:
            scene_dict[id(self)]["init_orientation"] = self._visualization_matrix[0]
        except:
            raise RuntimeError("Cannot generate visualization data " + \
                                "because numerical transformation " + \
                               "has not been performed, " + \
                                "Please call the numerical " + \
                               "transformation methods, " + \
                               "before generating visualization dict")

        return scene_dict

    def generate_simulation_dict(self):
        """
        Generates the simulation information for this visualization
        frame. It maps the simulation data information to the
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
            raise RuntimeError("Cannot generate visualization data "
                               "because numerical transformation "
                               "has not been performed, "
                               "Please call the numerical "
                               "transformation methods, "
                               "before generating visualization dict")

        return simulation_dict

    def _create_keyframetrack(self, times, dynamic_values, constant_values,
                              constant_map=None):
        """Sets attributes with a Mesh and KeyframeTrack for animating this
        visualization frame.

        Parameters
        ==========
        times : ndarray, shape(n,)
            Array of monotonically increasing or decreasing values of time.
        dynamics_values : ndarray, shape(n, m)
            Array of state values for each time.
        constant_values : array_like, shape(p,)
            Array of values for the constants.
        constant_map : dictionary
            A key value pair mapping from SymPy symbols to floating point
            values.

        Returns
        =======
        track : VectorKeyframeTrack
            PyThreeJS animation track.

        """
        # TODO : Passing in constant_values and constant_map is redundant,
        # right?
        if p3js is None:
            raise ImportError('pythreejs must be installed.')

        # NOTE : WebGL doesn't like 64bit so convert to 32 bit.
        times = np.asarray(times, dtype=np.float32)

        self._mesh = self.shape._p3js_mesh(constant_map=constant_map)

        # NOTE : This is required to set the transform matrix directly.
        self._mesh.matrixAutoUpdate = False

        matrices = self.evaluate_transformation_matrix(dynamic_values,
                                                       constant_values)

        # Set initial configuration.
        self._mesh.matrix = matrices[0]

        # TODO : If the user does not name their shapes, then there will be
        # KeyFrameTracks with duplicate names. Need a better fix for this, but
        # I at least warn the user if they didn't change the name at all.
        if self._mesh.name == 'unnamed':
            msg = ("The shape provided to this visualization frame must have a "
                   "unique name other thane 'unnamed'. Make sure all shapes in "
                   "the scene have unique names.")
            raise ValueError(msg)

        name = "scene/{}.matrix".format(self._mesh.name)

        track = p3js.VectorKeyframeTrack(name=name, times=times,
                                         values=matrices)

        self._track = track
