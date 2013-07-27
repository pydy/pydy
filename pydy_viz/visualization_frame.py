from sympy.physics.mechanics import *
from shapes import *
#TODO Docstrings

class VisualizationFrame(object):
    """
    A VisualizationFrame is an object used to draw a particular shape
    or set of shapes. 
    It is an implementation of Scene-Graph. 
   
    """
    def __init__(self, *args,**kwargs):
        #So, We use this approach for smart arguments..
        i=0
        #If args[i] matches a condition, it is put up and i is 
        #incremented ..
        
        #With Great power comes great responsibility!!
        
        #If nothing is supplied ...
        if not args and not kwargs:
            print '''WARNING:Initializing VisualizationFrame 
                              without any ReferenceFrame, Point.
                              For Visualizations to work, 
                              You must supply a ReferenceFrame, and a
                              Point'''
            self._name = 'UnNamed'
        
        #If name is supplied ...
        if i < len(args) and isinstance(args[i],str):
            self._name = args[0]
            i+=1
        #otherwise UnNamed ...    
        else:
            self._name = 'UnNamed'
            
        #RigidBody is supplied ...
        if i < len(args) and isinstance(args[i], RigidBody):
            self._reference_frame = args[i].get_frame()
            self._origin = args[i].get_masscenter()
            i+=1
            
        #Point,referenceFrame is supplied ...     
        if i+1<len(args) and isinstance(args[i], Point) and \
               isinstance(args[i+1], ReferenceFrame):
                   
            self._origin = args[i]
            self._reference_frame = args[i+1]
            i+=2    
        
        #ReferenceFrame/Point is supplied ...
        if i+1<len(args) and isinstance(args[i], ReferenceFrame) and \
               isinstance(args[i+1], Point):
                   
            self._reference_frame = args[i]
            self._origin = args[i+1]
            i+=2    
            
        if i<len(args) and isinstance(args[i], Shape):
            self._shape = args[i]
            i+=1    
        
        #User might want to use keyword args for assigning stuff ...
        if kwargs.has_key('name'):
            self._name = kwargs['name']
            
        if kwargs.has_key('reference_frame'):
            self._reference_frame = kwargs['reference_frame']
            
        if kwargs.has_key('origin'):
            self._origin = args['origin']        
            i+=1        
            
        if kwargs.has_key('shape'):
            self._shape = kwargs['shape']

    def transform(self, reference_frame, point):
        _rotation_matrix = self._reference_frame.dcm(reference_frame)

        self._transform[0:3, 0:3] = _rotation_matrix[0:3, 0:3]

        _point_vector = self._origin.pos_from(point).express(reference_frame)

        self._transform[3, 0] = _point_vector.dot(reference_frame.x)
        self._transform[3, 1] = _point_vector.dot(reference_frame.y)
        self._transform[3, 2] = _point_vector.dot(reference_frame.z)

        return self._transform
        
    def generate_numeric_transform(self, dynamic, parameters):
        """Returns a function which returns a transformation matrix given
        the symbolic states and the symbolic system parameters."""

        dummy_symbols = [Dummy() for i in dynamic]
        dummy_dict = dict(zip(dynamic, dummy_symbols))
        transform = self._transform.subs(dummy_dict)

        self.numeric_transform = lambdify(dummy_symbols + parameters,
                                          transform, modules="numpy")
   
   def evaluate_numeric_transform(self, states, parameters):
       #What if there is no such states.shape attrib?
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

        if len(states.shape) > 1:
            n = states.shape[0]
            new = zeros((n, 4, 4))
            for i, time_instance in enumerate(states):
                args = hstack((time_instance, parameters))
                new[i, :, :] = self.numeric_transform(*args)

        else:
            args = hstack((states, parameters))
            new = self.numeric_transform(*args)

        self.simulation_matrix = new
   def generate_simulation_dict(self):
        self._data = {}
        self._data['name'] = self._name
        self._data['shape'] = {}
        self._data['shape'] = self._shape.generate_data()  # hidden method
        self._data['simulation_matrix'] = self.simulation_matrix.tolist()
        return self._data 
