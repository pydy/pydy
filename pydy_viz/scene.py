#TODO imports

class Scene():
    """
    Scene class holds all the data required for the visualizations/
    animation of a system.
    
    It has methods for inputting the numerical data from the numerical 
    integrations of Equations of Motions and convert them to JSON 
    values, which can be then parsed by Javascripts(webgls).
    
    A scene object takes a ReferenceFrame, and a Point as required
    arguments. The reference_frame and point act as the inertial
    frame and origin with respect to which all objects are oriented
    and rendered in the visualizations
    
    A scene needs to be supplied with visualization_frames, Cameras, 
    and Light objects, as optional arguments.
    A scene can also be supplied with the height and width of the 
    browser window where visualization would be displayed.
    Default is 800 * 800.
    
    
    """
    
    def __init_(self, *args, width=800, height=800, \
                        visualization_frames=[], cameras=[], lights=[]):
        """
        Initializes a Scene instance.
        It requires a reference frame and a point to be initialized.
        
        Parameters:
        ===========
        
        name : str
        Name of Scene object.
        
        reference_frame : ReferenceFrame
        All the transformations would be carried out with respect
        to this reference frame. 
        
        origin : Point
        
        
        """
        i=0
        if isinstance(args[0], str):
            self._name = args[0]
            i+=1
        else:
            self._name = 'UnNamed'
        #DO we need isintance checks here?
        #To check valid rframe, point are supplied?
        
        self._reference_frame = args[i]
        self._origin = args[i+1]        
            
            
            
        #Setting cameras first ...
        #If not supplied ... 
        if not cameras:
            self.cameras = [PerspectiveCamera(self.reference_frame, \
                                                           self.point)]
        else:
            self.cameras = cameras
        if not lights:
            #self.lights = [] TODO Lights
            pass
        else:
            self.lights = lights    
        
        #Do we need additional methods for adding_vframes?
        #Or just use lists append and remove methods?
        
        self.visualization_frames = visualization_frames
        self.width = width
        self.height = height
        
        @property
        def name(self):
            """
            Name of Scene.
            """
            return self._name
            
        @name.setter
        def name(self, new_name):
            if not isinstance(new_name, str):
                raise TypeError('Name should be a valid str.')
            else:
                self._name = new_name
        
           @property
    def origin(self):
        """
        origin of the Scene.
        
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
        reference_frame of the Scene.
        
        """
        return self._reference_frame
        
    @reference_frame.setter
    def reference_frame(self, new_reference_frame):
        if not isinstance(new_reference_frame, ReferenceFrame):
            raise TypeError('''reference_frame should be a valid
                                ReferenceFrame object.''')
        else:
            self._reference_frame = new_reference_frame             
    
    def generate_json(self, state_sym, par_sym, states, parameters, \
                                              save_to='data.json'):
        """
        TODO: Docstring
        """
        out_file = save_to
        self._scene_data = {}
        self._scene_data['name'] = self._name
        self._scene_data['height'] = self._height
        self._scene_data['width'] = self._width
        self._scene_data['frames'] = []

        for frame in self.visualization_frames:

            frame.transform(self._reference_frame, self._origin)
            frame.generate_numeric_transform(state_sym, par_sym)
            frame.evaluate_numeric_transform(states, parameters)
            self._scene_data['frames'].append(frame.generate_simulation_dict())

        out_file = open()
        #return self._scene_data
        

    def display(self, json=None):
        #TODO
        
