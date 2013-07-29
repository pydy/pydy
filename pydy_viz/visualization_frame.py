from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, \
                                  Point, RigidBody, Particle, inertia

from shapes import Shape
import numpy as np
from sympy.matrices.expressions import Identity
from sympy import Dummy, lambdify

class VisualizationFrame(object):
    """
    A VisualizationFrame is an object used to draw a particular shape.
    It provides a ReferenceFrame, and an origin to the
    particular shape attached to it, and is hence useful for the 
    visualization and animation of the Shape.
    
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
    
    
    Parameters
    ==========
    name : str
        Name assigned to VisualizationFrame, default is UnNamed
    shape : Shape
        Shape to be attached to the VisualizationFrame
    reference_frame : ReferenceFrame
        reference_frame with respect to which all orientations of the 
        shape takes place, during visualizations/animations.
    origin : Point
        point with respect to which all the translations of the shape 
        takes place, during visualizations/animations.
       
    """
    def __init__(self, *args):
        """
        Initialises a VisualizationFrame object.
        To initialize a visualization frame, we need to supply
        a name(optional), a reference frame, a point, a shape. 
        
        Examples
        ========
        >>> from pydy_viz import VisualizationFrame, Shape
        >>> from sympy.physics.mechanics import \
                               ReferenceFrame, Point, RigidBody, \
                                Particle, inertia
        >>> from sympy import symbols                               
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')
        >>> shape = Shape()
        >>> #initializing with reference frame, point
        >>> frame1 = VisualizationFrame('frame1', I, O, shape)
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O)
        >>> # Initializing with a rigidbody ..
        >>> frame2 = VisualizationFrame('frame2', rbody, shape)
        >>> Pa = Particle('Pa', O, mass)
        >>> #initializing with Particle, reference_frame ...
        >>> frame3 = VisualizationFrame('frame3', I, Pa, shape)
        >>> #These can be changed later too ..
        >>> frame1.name = 'frame1_'
        >>> frame1.name
        'frame1_'
        >>> frame1.reference_frame = I
        >>> frame1.reference_frame
        I
        >>> frame1.shape = shape
        >>> frame1.shape
        shape
        """
        #Last arg should be a Shape ..
        if isinstance(args[-1], Shape):
            self._shape = args[-1]
        else:
            raise TypeError('''Please provide a valid shape object''')    
        i = 0
        #If first arg is not str, name the visualization frame 'UnNamed'    
        if isinstance(args[i], str):
            self._name = args[i]
            i += 1
        else:
            self._name = 'UnNamed'        
        
        try:
            self._reference_frame = args[i].get_frame()
            self._origin = args[i].get_masscenter()
            
        except AttributeError:
            #It is not a rigidbody, hence this arg should be a 
            #reference frame
            self._reference_frame = args[i]
            i += 1
            
            #Now next arg can either be a Particle or point
            try:
                self._origin = args[i].get_point()
            except AttributeError:
                self._origin = args[i]
                
        #basic things required, transform matrix and child frames  
        self._transform = Identity(4).as_mutable()
        self.child_frames = []
        
    #setting attributes ..
    def __str__(self):
        return 'VisualizationFrame ' + self._name
        
    def __repr__(self):
        return 'VisualizationFrame'    
        
    @property
    def name(self):
        """
        name attribute of the Visualization Frame
        """
        return self._name
    @name.setter
    def name(self,new_name):
        if not isinstance(new_name, str):
            raise TypeError('''Name should be a str object''')
        else:
            self._name = new_name
    
    @property
    def origin(self):
        """
        origin attribute of the VisualizationFrame, 
        with respect to which all translational transformations
        take place.
        """
        return self._origin
    @origin.setter
    def origin(self, new_origin):
        if not isinstance(new_origin, Point):
            raise TypeError('''origin should be a valid Point Object''')
        else:
            self._origin = new_origin        
            
    @property
    def reference_frame(self):
        """
        reference_frame attribute of the VisualizationFrame, 
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
        shape attribute of the VisualizationFrame.
        A shape attached to the visualization frame.
        NOTE: Only one shape can be attached to a visualization frame.
        """
        return self._shape
    @shape.setter
    def shape(self,new_shape):
        if not isinstance(new_shape, Shape):
            raise TypeError('''shape should be a valid Shape object.''')    
        else:
            self._shape = new_shape
               
    def transform(self, reference_frame, point):
        _rotation_matrix = self.reference_frame.dcm(reference_frame)

        self._transform[0:3, 0:3] = _rotation_matrix[0:3, 0:3]

        _point_vector = self.origin.pos_from(point).express(reference_frame)

        self._transform[3, 0] = _point_vector.dot(reference_frame.x)
        self._transform[3, 1] = _point_vector.dot(reference_frame.y)
        self._transform[3, 2] = _point_vector.dot(reference_frame.z)

        return self._transform
        
    def generate_numeric_transform(self, dynamic, parameters):
        """Returns a function which returns a transformation matrix given
        the symbolic states and the symbolic system parameters.
        
        Parameters
        ==========
        dynamic : list of all the dynamic symbols used in defining the 
                  mechanics objects.
        parameters : list of all symbols used in defining the 
                     mechanics objects
    
        """

        dummy_symbols = [Dummy() for i in dynamic]
        dummy_dict = dict(zip(dynamic, dummy_symbols))
        transform = self._transform.subs(dummy_dict)

        self.numeric_transform = lambdify(dummy_symbols + parameters,
                                          transform, modules="numpy")
   
    def evaluate_numeric_transform(self, states, parameters):
        """Returns the numerical transformation matrices for each time step.

        Parameters
        ----------
        states : array_like, shape(m,) or shape(n, m)
            The m states for each n time step.
        parameters : array_like, shape(p,)
            The p constant parameters of the system.

        Returns
        -------
        transform_matrix : numpy.array, shape(n, 4, 4)
            A 4 x 4 transformation matrix for each time step.

        """
        #If states is instance of numpy array, well and good.
        #else convert it to one:
        if not isinstance(states,np.ndarray):
            states = np.array(states)
        if len(states.shape) > 1:
            n = states.shape[0]
            new = zeros((n, 4, 4))
            for i, time_instance in enumerate(states):
                args = np.hstack((time_instance, parameters))
                new[i, :, :] = self.numeric_transform(*args)

        else:
            args = np.hstack((states, parameters))
            new = self.numeric_transform(*args)

        self.simulation_matrix = new
        return self.simulation_matrix
        
    def add_child_frames(self, *args):
        """
        Used for nesting of visualization frames.
        Helpful for implementing Scene graph.
        
        In this way we can create complex shapes.
        We break complex shapes into simple shapes from the Shape
        subclasses, then we supply a visualization frame for them,
        and then we add them to a parent frame, whose reference frame
        and origin are those in which the shape is created on.
        
        Parameters
        ----------
        visualization_frame: one or more VisualizationFrame objects.

        """
        for _frame in args:
            if not isinstance(_frame, VisualizationFrame):
                raise TypeError('''frame is not an instance 
                                   of VisualizationFrame:''' + _frame)    
            else:
                self.child_frames.append(_frame)      
                                 
                    
    def generate_simulation_dict(self):
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
        self._data = {}
        self._data['name'] = self.name
        self._data['children'] = []
        
        #if child_frames list isnt empty
        if len(self.child_frames > 0):
            for _child in self.child_frames:
                self._data['children'].append( \
                                     _child.generate_simulation_dict())
             
        self._data['shape'] = self.shape.generate_data()  
            
        if not self.simulation_matrix:
            raise Error('''Please call the numerical 
                            transformation methods,
                           before generating simulation dict ''')
        else:                   
            self._data['simulation_matrix'] = self.simulation_matrix.tolist()
            
        return self._data 

