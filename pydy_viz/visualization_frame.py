from sympy.physics.mechanics import *
from shapes import *
import numpy as np
from sympy.matrices.expressions import Identity
from sympy import Dummy, lambdify

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
            self._name = args[i]
            i+=1
        #otherwise UnNamed ...    
        else:
            self._name = 'UnNamed'
            
        #RigidBody is supplied ...
        if i < len(args) and isinstance(args[i], RigidBody):
            self._reference_frame = args[i].get_frame()
            self._origin = args[i].get_masscenter()
            i+=1
            
        #Point(or)Particle,referenceFrame is supplied ...     
        if i+1<len(args):
            if isinstance(args[i], Point) and \
               isinstance(args[i+1], ReferenceFrame):
                   
                self._origin = args[i]
                self._reference_frame = args[i+1]
                i+=2
                
            elif isinstance(args[i], Particle) and \
               isinstance(args[i+1], ReferenceFrame):
                   
                self._origin = args[i].get_point()
                self._reference_frame = args[i+1]
                i+=2        
        
        #ReferenceFrame,Point(or Particle) is supplied ...
        if i+1<len(args): 
            if isinstance(args[i], ReferenceFrame) and \
               isinstance(args[i+1], Point):
                       
                self._reference_frame = args[i]
                self._origin = args[i+1]
                i+=2    
            elif isinstance(args[i], ReferenceFrame) and \
                 isinstance(args[i+1], Particle):
                   
                self._reference_frame = args[i]
                self._origin = args[i+1].get_point()
                i+=2    
            
        if i<len(args) and isinstance(args[i], Shape):
            self._shape = args[i]
            i+=1    
        
        #User might want to use keyword args for assigning stuff ...
        if kwargs.has_key('name'):
            if not isinstance(kwargs['name'], str):
                raise TypeError('''Name should be a str object''')
            else:
                self._name = kwargs['name']
            
        if kwargs.has_key('reference_frame'):
            if not isinstance(kwargs['reference_frame'], \
                                             ReferenceFrame):
                raise TypeError('''reference_frame should be a valid
                                         ReferenceFrame object''')
            else:
                self._reference_frame = kwargs['reference_frame']
            
        if kwargs.has_key('origin'):
            if not isinstance(kwargs['origin'], \
                                             Point):
                raise TypeError('''origin should be a valid
                                         Point object''')
            else:
                self._origin = kwargs['origin']
            
        if kwargs.has_key('shape'):
            if not isinstance(kwargs['shape'], \
                                             Shape):
                raise TypeError('''shape should be a valid
                                         Shape object''')
            else:
                self._shape = kwargs['shape']
                
        if kwargs.has_key('particle'):
            if not isinstance(kwargs['particle'], \
                                             Particle):
                raise TypeError('''particle should be a valid
                                         Particle object''')
            else:
                self._origin = kwargs['particle'].get_point        
                
        if kwargs.has_key('rigidbody'):
            if not isinstance(kwargs['rigidbody'], \
                                             RigidBody):
                raise TypeError('''rigidbody should be a valid RigidBody''')
                
            else:
                self._reference_frame = kwargs['rigidbody'].get_frame()
                self._origin = \
                                   kwargs['rigidbody'].get_masscenter()
                                   
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
        return self._name
    @name.setter
    def name(self,new_name):
        if not isinstance(new_name, str):
            raise TypeError('''Name should be a str object''')
        else:
            self._name = new_name
    
    @property
    def origin(self):
        return self._origin
    @origin.setter
    def origin(self, new_origin):
        if not isinstance(new_origin, Point):
            raise TypeError('''origin should be a valid Point Object''')
        else:
            self._origin = new_origin        
            
    @property
    def reference_frame(self):
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
        return self._shape
    @shape.setter
    def shape(self,new_shape):
        if not isinstance(new_shape, Shape):
            raise TypeError('''shape should be a valid Shape object.''')    
        else:
            self._shape = new_shape
               
    def transform(self, reference_frame, point):
        if not self.reference_frame:
            raise Error('''Please supply a reference_frame to 
                            the VisualizationFrame instance.''')
        if not self.origin:
            raise Error('''Please supply an origin to 
                            the VisualizationFrame instance.''')                    
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
        for _frame in args:
            if not isinstance(_frame, VisualizationFrame):
                raise TypeError('''frame is not an instance 
                                   of VisualizationFrame:''' + _frame)    
            else:
                self.child_frames.append(_frame)      
                                 
    def remove_child_frames(self, *args):
        for _frame in args:
            if not isinstance(_frame, VisualizationFrame):
                    raise TypeError('''frame is not an instance 
                                       of VisualizationFrame:''' + _frame)    
            else:
                self.child_frames.remove(_frame)          
                    
    def generate_simulation_dict(self):
        self._data = {}
        self._data['name'] = self.name
        self._data['children'] = []
        
        #if child_frames list isnt empty
        if len(self.child_frames > 0):
            for _child in self.child_frames:
                self._data['children'].append( \
                                     _child.generate_simulation_dict())
        if not self.shape:
            raise Error('''Please assign a shape to the 
                          VisualizationFrame before calling 
                          generate_simulation_dict()''')
        else:                  
            self._data['shape'] = self.shape.generate_data()  
            
        if not self.simulation_matrix:
            raise Error('''Please call the numerical 
                            transformation methods,
                           before generating simulation dict ''')
        else:                   
            self._data['simulation_matrix'] = self.simulation_matrix.tolist()
            
        return self._data 

