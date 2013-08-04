from sympy.physics.mechanics import 
import json

try:
    import IPython
except:
    pass
        
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
        
        self._reference_frame = args[i]
        self._origin = args[i+1]        
            
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
        Origin of the Scene.
        
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
        generate_json() method generates a dictionary, which is dumped
        as JSON in the file name given by save_to argument.
        Default filename is data.json.
        The JSON file contain following keys:
        1) Width of the scene.
        2) Height of the scene.
        3) name of the scene.
        4) frames in the scene, which contains sub-dictionaries 
           of all the visualization frames information.
        

        Parameters
        ==========
        state_sym : Sympifyable list or tuple
            This contains all the dynamic symbols or state variables
            which are required for solving the transformation matrices
            of all the frames of the scene.
       
        par_sym : Sympifyable list or tuple
            This contains all the symbols for the parameters which are
            used for defining various objects in the system.
     
        states : list or tuple
            initial states of the system. The list or tuple 
            should be respective to the state_sym.

        parameters : list or tuple
            values of the parameters. The list or tuple 
            should be respective to the par_sym.

        save_to : str
            path to the file where to write the generated data JSON.
            the path should be chosen such as to have the write 
            permissions to the user.

        Examples
        ========
        #TODO : Write complete example for initializing a vframe and stuff, 
               or directly show this method?              
        """
        self.saved_json_file = save_to
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

     
        outfile = open(self.saved_json_file)
        outfile.write(json.dumps(self._scene_data, indent=4, separators=(',', ': '))) 
        outfile.close()

    def _display_from_ipython(json_data=None):
        #Basic script requirements. ..
            

    def _display_from_interpreter(json_data=None):     
        #Get html file ..
        _path_to_html = pydy_viz.__file__[:-12] + 'index.html'
        content = open(_path_to_html)
 
        #If json_data is not provided, use self.saved_json_file
        #Replace json object name in content to self.saved_json_file 
        _json_replacement = json_data or self.saved_json_file
        content = 
        
        out_file = open('index.html','w')
        outfile.write(content)
        
        start_server()
 


    def display(self, json_data=None):
        """
        display method can be used in two ways.
        When called from IPython notebook, it shows the visualization
        in the form of output cell in the IPython notebook.
        If it is called from python interpreter or 
        IPython interpreter(not notebook), It generates an html file,
        in the current directory, which can be opened in the webgl 
        compliant browser for viewing the visualizations.
        
        This method can also be used to load any json file, irrespective
        of whether it was created in the same session,
        
        Parameters
        ==========
        json_data : str
            path to the json file which is to be visualized.
            (optional).
        
        Examples
        =======

        """        
     
        try:
        config = get_ipython().config
            if config['KernelApp']['parent_appname'] == 'ipython-notebook':
                #Launched from IPython browser
                pass
        else:
            #Launched from interpreter        


        
