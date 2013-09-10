__all__ = ['VisualizationFrame']

from sympy.physics.mechanics import Point, ReferenceFrame
from shapes import Shape
import numpy as np
from sympy.matrices.expressions import Identity
from sympy import Dummy, lambdify


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
        >>> from pydy_viz import VisualizationFrame, Sphere
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
            raise TypeError('''Please provide a valid shape object''')
        i = 0
        #If first arg is not str, name the visualization frame 'unnamed'
        if isinstance(args[i], str):
            self._name = args[i]
            i += 1
        else:
            self._name = 'unnamed'

        try:
            self._reference_frame = args[i].get_frame()
            self._origin = args[i].get_masscenter()

        except AttributeError:
            #It is not a rigidbody, hence this arg should be a
            #reference frame
            try:
                dcm = args[i].dcm()
                self._reference_frame = args[i]
                i += 1
            except AttributeError:
                raise TypeError(''' A ReferenceFrame is to be supplied
                                   before a Particle/Point. ''')    

            #Now next arg can either be a Particle or point
            try:
                self._origin = args[i].get_point()
                
            except AttributeError:
                self._origin = args[i]

    #setting attributes ..
    def __str__(self):
        return 'VisualizationFrame ' + self._name

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
        """
        Generates the symbolic Transformation matrix,
        with respect to the reference_frame, point in the argument.

        Parameters
        ==========

        reference_frame : ReferenceFrame
        A reference_frame with respect to which transformation matrix
        is generated.
        point : Point
        A point with respect to which transformation matrix
        is generated.

        Returns
        =======
        A SymPy 4by4 matrix, containing symbolic variables for
        transformation.

        """
        _rotation_matrix = self.reference_frame.dcm(reference_frame)
        self._transform = Identity(4).as_mutable()
        self._transform[0:3, 0:3] = _rotation_matrix[0:3, 0:3]

        print self.origin.pos_from(point)
        print type(self.origin.pos_from(point))
        _point_vector = self.origin.pos_from(point).express(reference_frame)

        self._transform[3, 0] = _point_vector.dot(reference_frame.x)
        self._transform[3, 1] = _point_vector.dot(reference_frame.y)
        self._transform[3, 2] = _point_vector.dot(reference_frame.z)
        return self._transform

    def generate_numeric_transform_function(self, dynamic_variables, constant_variables):
        """Returns a function which returns a transformation matrix given
        the symbolic states and the symbolic system parameters.


        Parameters
        ==========
        dynamic_variables : list of all the dynamic symbols used in defining the
                  mechanics objects.
        constant_variables : list of all symbols used in defining the
                     mechanics objects

        Returns
        =======
        A Lambda function which returns a transformation matrix,
        given symbolic states, and symbolic system parameters

        """

        dummy_symbols = [Dummy() for i in dynamic_variables]
        dummy_dict = dict(zip(dynamic_variables, dummy_symbols))
        transform = self._transform.subs(dummy_dict)

        self._numeric_transform = lambdify(dummy_symbols + constant_variables,
                                          transform, modules="numpy")

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
        transform_matrix : numpy.array, shape(n, 4, 4)
            A 4 x 4 transformation matrix for each time step.

        """
        #If states is instance of numpy array, well and good.
        #else convert it to one:

        states = np.array(dynamic_values)
        if len(states.shape) > 1:
            n = states.shape[0]
            new = np.zeros((n, 4, 4))
            for i, time_instance in enumerate(states):
                args = np.hstack((time_instance, constant_values))
                new[i, :, :] = self._numeric_transform(*args)

        else:
            n = 1
            args = np.hstack((states, constant_values))
            new = self._numeric_transform(*args)

        print 'N:', n
        print 'new:', new
        self._visualization_matrix = new.reshape(n, 16)
        print 'after reshape:', self._visualization_matrix
        return self._visualization_matrix

    def generate_visualization_dict(self):
        """
        Returns a dictionary of all the info required
        for the visualization of this frame, alongwith child.

        Before calling this method, all the transformation matrix
        generation methods should be called, or it will give an error.

        Returns
        ======

        a dictionary containing following keys:

        name : name of the VisualizationFrame
        children : simulation dictionary of child frames
        shape : shape info of the attached shape,
        like dimensions, color etc.It is generated from generator method
        of Shape class.
        simulation_matrix : a N*4*4 matrix, converted to list, for
        passing to Javascript for animation purposes, where N is the
        number of timesteps for animations.
        """
        _data = {}
        _data['name'] = self.name
        _data['shape'] = self.shape.generate_dict()
        try:
            _data['simulation_matrix'] = \
                                     self._visualization_matrix.tolist()
        except:
            raise RuntimeError('''Please call the numerical
                            transformation methods,
                           before generating visualization dict ''')


        return _data
