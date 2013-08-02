
from visualization_frame import VisualizationFrame

class PerspectiveCamera(VisualizationFrame):
    """
    Creates a Perspective Camera for visualization. 
    The camera is inherited from VisualizationFrame,
    
    It can be attached to dynamics objects, hence we can 
    get a moving camera. All the transformation matrix generation
    methods are applicable to a Perspective Camera.
    Like VisualizationFrame,
    It can also be initialized using:
    1)Rigidbody
    2)ReferenceFrame, Point
    3)ReferenceFrame, Particle
    Either one of these must be supplied during initialization
    
    Unlike VisualizationFrame, It doesnt require a Shape argument.
    
    Parameters:
    ===========
    
    name : str
    a name for the PerspectiveCamera(optional). Default is 'UnNamed'
    
    fov : int or float
    Field Of View, It determines the angle between the top and bottom 
    of the viewable area(in degrees). Default is 45 (degrees)
    
    near : int or float
    The distance of near plane of the PerspectiveCamera.
    All objects closer to this distance are not displayed.
    
    far : int or float
    The distance of far plane of the PerspectiveCamera
    All objects farther than this distance are not displayed.
    
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialises a PerspectiveCamera object.
        To initialize a visualization frame, we need to supply
        a name(optional), a reference frame, a point, 
        field of view(fov) (optional), near plane distance(optional)
        and far plane distance(optional).
        
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
        >>> camera1 = PerspectiveCamera('frame1', I, O)
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> # Initializing with a rigidbody ..
        >>> camera2 = PerspectiveCamera('frame2', rbody)
        >>> Pa = Particle('Pa', O, mass)
        >>> #initializing with Particle, reference_frame ...
        >>> camera3 = PerspectiveCamera('frame3', I, Pa)
        """
        try:
            self._fov = kwargs['fov']
        except KeyError:
            self._fov = 45    
      
        try:
            self._near = kwargs['near']
        except KeyError:
            self._near = 1
            
        try:
            self._far = kwargs['far']
        except KeyError:
            self._far = 1000    
                
        
        #Now we use same approach as in VisualizationFrame
        #for setting reference_frame and origin 
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
                
        #basic thing required, transform matrix
        self._transform = Identity(4).as_mutable()

    def __str__(self):
        return 'PerspectiveCamera: ' + self._name
    
    def __repr__(self):
        return 'PerspectiveCamera'             
        
    @property
    def fov(self):
        """
        attribute for Field Of view of a PerspectiveCamera
        Default is 45 degrees
        """
        return self._fov
    @fov.setter
    def fov(self, new_fov):
        if not isinstance(new_fov, (int, str)):
            raise TypeError(''' fov should be supplied in 
                                         int or float ''')    
        else:
            self._fov = new_fov

    @property
    def near(self):
        """
        attribute for Near Plane distance of a PerspectiveCamera
        Default is 1
        """
        return self._near
    @near.setter
    def near(self, new_near):
        if not isinstance(new_near, (int, str)):
            raise TypeError(''' near should be supplied in 
                                         int or float ''')    
        else:
            self._near = new_near

    @property
    def far(self):
        """
        attribute for Far Plane distance of a PerspectiveCamera
        Default is 1000
        """
        return self._far
    @far.setter
    def far(self, new_far):
        if not isinstance(new_far, (int, str)):
            raise TypeError(''' far should be supplied in 
                                         int or float ''')    
        else:
            self._far = new_far

    def generate_simulation_dict(self):
        """
        Returns a dictionary of all the info required
        for the visualization of this Camera
        
        Before calling this method, all the transformation matrix
        generation methods should be called, or it will give an error.
        
        Returns
        ======
        
        a dictionary containing following keys:
        
        name : name of the PerspectiveCamera
        fov : Field Of View of the PerspectiveCamera
        simulation_matrix : a N*4*4 matrix, converted to list, for
        passing to Javascript for animation purposes, where N is the 
        number of timesteps for animations.
        
        
        """
        self._data = {}
        self._data['name'] = self.name
        self._data['type'] = self.__repr__()
        self._data['fov'] = self.fov
        self._data['near'] = self.near
        self._data['far'] = self.far
        
        if not self.simulation_matrix:
            #Not sure which error to call here.
            raise RuntimeError('''Please call the numerical 
                            transformation methods,
                           before generating simulation dict ''')
        else:                   
            self._data['simulation_matrix'] = self.simulation_matrix.tolist()
            
        return self._data 

class OrthoGraphicCamera(VisualizationFrame):
    """
    Creates a OrthoGraphic Camera for visualization. 
    The camera is inherited from VisualizationFrame,
    
    It can be attached to dynamics objects, hence we can 
    get a moving camera. All the transformation matrix generation
    methods are applicable to a Perspective Camera.
    Like VisualizationFrame,
    It can also be initialized using:
    1)Rigidbody
    2)ReferenceFrame, Point
    3)ReferenceFrame, Particle
    Either one of these must be supplied during initialization
    
    Unlike VisualizationFrame, It doesnt require a Shape argument.
    
    Parameters:
    ===========
    
    name : str
    a name for the PerspectiveCamera(optional). Default is 'UnNamed'
    
    near : int or float
    The distance of near plane of the PerspectiveCamera.
    All objects closer to this distance are not displayed.
    
    far : int or float
    The distance of far plane of the PerspectiveCamera
    All objects farther than this distance are not displayed.
    
    """
    
    def __init__(self, *args, **kwargs):
        """
        Initialises an OrthoGraphicCamera object.
        To initialize a visualization frame, we need to supply
        a name(optional), a reference frame, a point, 
        near plane distance(optional)
        and far plane distance(optional).
        
        Examples
        ========
        >>> from pydy_viz import OrthoGraphicCamera
        >>> from sympy.physics.mechanics import \
                               ReferenceFrame, Point, RigidBody, \
                                Particle, inertia
        >>> from sympy import symbols                               
        >>> I = ReferenceFrame('I')
        >>> O = Point('O')
        >>> shape = Shape()
        >>> #initializing with reference frame, point
        >>> camera1 = OrthoGraphicCamera('frame1', I, O)
        >>> Ixx, Iyy, Izz, mass = symbols('Ixx Iyy Izz mass')
        >>> i = inertia(I, Ixx, Iyy, Izz)
        >>> rbody = RigidBody('rbody', O, I, mass, (inertia, O))
        >>> # Initializing with a rigidbody ..
        >>> camera2 = OrthoGraphicCamera('frame2', rbody)
        >>> Pa = Particle('Pa', O, mass)
        >>> #initializing with Particle, reference_frame ...
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
                
        
        #Now we use same approach as in VisualizationFrame
        #for setting reference_frame and origin 
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
                
        #basic thing required, transform matrix
        self._transform = Identity(4).as_mutable()
         
    def __str__(self):
        return 'OrthoGraphicCamera: ' + self._name
    
    def __repr__(self):
        return 'OrthoGraphicCamera'             
    
    @property
    def near(self):
        """
        attribute for Near Plane distance of an OrthoGraphicCamera
        Default is 1
        """
        return self._near
    @near.setter
    def near(self, new_near):
        if not isinstance(new_near, (int, str)):
            raise TypeError(''' near should be supplied in 
                                         int or float ''')    
        else:
            self._near = new_near

    @property
    def far(self):
        """
        attribute for Far Plane distance of an OrthoGraphicCamera
        Default is 1000
        """
        return self._far
    @far.setter
    def far(self, new_far):
        if not isinstance(new_far, (int, str)):
            raise TypeError(''' far should be supplied in 
                                         int or float ''')    
        else:
            self._far = new_far

    def generate_simulation_dict(self):
        """
        Returns a dictionary of all the info required
        for the visualization of this Camera
        
        Before calling this method, all the transformation matrix
        generation methods should be called, or it will give an error.
        
        Returns
        ======
        
        a dictionary containing following keys:
        
        name : name of the OrthoGraphicCamera
        simulation_matrix : a N*4*4 matrix, converted to list, for
        passing to Javascript for animation purposes, where N is the 
        number of timesteps for animations.
        
        
        """
        self._data = {}
        self._data['name'] = self.name
        self._data['type'] = self.__repr__()
        self._data['near'] = self.near
        self._data['far'] = self.far
        
        if not self.simulation_matrix:
            #Not sure which error to call here.
            raise RuntimeError('''Please call the numerical 
                            transformation methods,
                           before generating simulation dict ''')
        else:                   
            self._data['simulation_matrix'] = self.simulation_matrix.tolist()
            
        return self._data 
