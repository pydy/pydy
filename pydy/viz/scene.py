#!/usr/bin/env python

# standard library
import os
import json
import distutils
import distutils.dir_util
import webbrowser
import datetime
# external
from sympy.physics.mechanics import ReferenceFrame, Point

# local
from .camera import PerspectiveCamera
from .server import Server
from .light import PointLight

__all__ = ['Scene']

try:
    from IPython.lib import backgroundjobs as bg
except ImportError:
    IPython = None


class Scene(object):
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
    def __init__(self, reference_frame, origin, *visualization_frames,
                 **kwargs):
        """Initializes a Scene instance.

        Parameters
        ==========
        reference_frame : ReferenceFrame
            All the transformations would be carried out with respect to
            this reference frame.
        origin : Point
            All the transformations would be carried out with respect to
            this point.
        visualization_frames : VisualizationFrame
            One or more visualization frames which are to be visualized in
            the scene.
        name : str, optional
            Name of Scene object.
        width : int or float, optional
            The width of the canvas used for visualizations. Default is
            800px.
        height : int or float, optional
            Height of the canvas used for visualizations. Default is 800px.
        camera : Camera, optional
            The camera with which to display the object. Default is
            PerspectiveCamera, with reference_frame and origin same as
            defined for this scene.
        """

        self._reference_frame = reference_frame
        self._origin = origin

        self.visualization_frames = list(visualization_frames)

        try:
            self._name = kwargs['name']
        except KeyError:
            self._name = 'unnamed'

        try:
            self._width = kwargs['width']
        except KeyError:
            self._width = 800

        try:
            self._height = kwargs['height']
        except KeyError:
            self._height = 800

        try:
            self.cameras = kwargs['cameras']
        except KeyError:
            self.cameras = [PerspectiveCamera('DefaultCamera',
                            self._reference_frame,
                            self._origin.locatenew(
                                'p_camera',
                                10*self._reference_frame.z))]

        try:
            self.lights = kwargs['lights']
        except KeyError:
            self.lights = [PointLight('DefaultLight',
                           self._reference_frame,
                           self._origin.locatenew(
                               'p_light',
                               10*self._reference_frame.z))]

    @property
    def name(self):
        """
        Returns Name of Scene.
        """
        return self._name

    @name.setter
    def name(self, new_name):
        """
        sets name of scene.
        """
        if not isinstance(new_name, str):
            raise TypeError('Name should be a valid str.')
        else:
            self._name = new_name

    @property
    def origin(self):
        """
        returns Origin of the Scene.
        """
        return self._origin

    @origin.setter
    def origin(self, new_origin):
        """
        sets origin of the scene
        """
        if not isinstance(new_origin, Point):
            raise TypeError('''origin should be a valid Point Object''')
        else:
            self._origin = new_origin

    @property
    def reference_frame(self):
        """
        returns reference_frame of the Scene.
        """
        return self._reference_frame

    @reference_frame.setter
    def reference_frame(self, new_reference_frame):
        """
        Sets reference frame for the scene.
        """
        if not isinstance(new_reference_frame, ReferenceFrame):
            raise TypeError('''reference_frame should be a valid
                                ReferenceFrame object.''')
        else:
            self._reference_frame = new_reference_frame

    def generate_visualization_json(self, dynamic_variables,
                                    constant_variables, dynamic_values,
                                    constant_values, outfile_prefix=None):
        """
        generate_visualization_json() method generates a json str, which is
        saved to file.

        Parameters
        ==========
        dynamic_variables : Sympifyable list or tuple
            This contains all the dynamic symbols or state variables
            which are required for solving the transformation matrices
            of all the frames of the scene.

        constant_variables : Sympifyable list or tuple
            This contains all the symbols for the parameters which are
            used for defining various objects in the system.

        dynamic_values : list or tuple
            initial states of the system. The list or tuple
            should be respective to the state_sym.

        constant_values : list or tuple
            values of the parameters. The list or tuple
            should be respective to the par_sym.

        outfile_prefix : str
            A prefix to be put while saving the scene_desc
            and simulation_data files. Files will be named
            as `outfile_prefix_scene_desc.json` and 
            `outfile_prefix_simulation_data.json`. If not specified
            a timestamp shall be used as the prefix.
            
        Returns
        =======

        The dictionary contains following keys:
        1) Width of the scene.
        2) Height of the scene.
        3) name of the scene.
        4) frames in the scene, which contains sub-dictionaries
           of all the visualization frames information.


        """
        if outfile_prefix is None:
            outfile_prefix = "_".join(str(datetime.datetime.now()).\
                                  split(".")[0].split(" "))
            
        constant_map = dict(zip(constant_variables, constant_values))
        self.scene_json_file = outfile_prefix + "_scene_desc.json"
        self.simulation_json_file = outfile_prefix + "_simulation_data.json"

        self._simulation_data_dict = self.generate_simulation_dict(dynamic_variables,
                                                           constant_variables,
                                                           dynamic_values,
                                                           constant_values)
        self._scene_data_dict = self.generate_scene_dict(constant_map=constant_map)
        
        scene_data_outfile = open(self.scene_json_file, 'w')
        scene_data_outfile.write(json.dumps(self._scene_data_dict, indent=4,
                                 separators=(',', ': ')))
        scene_data_outfile.close()

        simulation_data_outfile = open(self.simulation_json_file, 'w')
        simulation_data_outfile.write(json.dumps(self._simulation_data_dict, indent=4,
                                 separators=(',', ': ')))
        simulation_data_outfile.close()

    def generate_scene_dict(self, constant_map={}):
        """
        This method is used to create the dictionary compatible with 
        MotionView. This JSON file contains all the relevant information
        required by MotionView to draw the scene on the canvas.
        
        
        """

        self._scene_info = {}
        self._scene_info["source"] = "PyDy"
        self._scene_info["name"] = self._name
        self._scene_info["newtonian_frame"] = str(self._reference_frame)
        self._scene_info["workspaceSize"] = 0.2#This should be accomodated in scene
                                                #instead of width/height of scene
        self._scene_info["objects"] = []
        
        for frame in self.visualization_frames:
            _object_info = frame.generate_scene_dict(constant_map=constant_map)
            self._scene_info["objects"].append(_object_info)        
        
        return self._scene_info
            
    def generate_simulation_dict(self, dynamic_variables,
                                    constant_variables, dynamic_values,
                                    constant_values):
        """
        This method is used to create the JSON file compatible with
        MotionView. This JSON file consists of all the simulation data
        along with references to the objects, for allowing motion to the
        objects in the MotionView visualizer.
        
        """
        #TODO: Add code to include lights and Cameras as well

        self._simulation_info = {}

        for frame in self.visualization_frames:
            frame.generate_transformation_matrix(self._reference_frame,
                                                 self._origin)
            frame.generate_numeric_transform_function(dynamic_variables,
                                                      constant_variables)
            frame.evaluate_transformation_matrix(dynamic_values, 
                                                         constant_values)

            self._simulation_info.update(frame.generate_simulation_dict())

        return self._simulation_info    
        
    def _copy_static_dir(self):
        """
        Copies all the static files required in the
        visualization to the current working directory
        The files and sub directories are stored within
        a hidden directory named .pydy_viz in the current
        working directory.
        Working directory can be cleaned by calling _cleanup()
        method, which deletes the .pydy_viz directory.

        """
        dst = os.path.join(os.getcwd(), '.pydy_viz')
        os.mkdir(dst)
        src = os.path.join(os.path.dirname(__file__), 'static')
        distutils.dir_util.copy_tree(src, dst)

    def _cleanup(self):
        distutils.dir_util.remove_tree(os.path.join(os.getcwd(), '.pydy_viz'))

    def _display_from_interpreter(self):
        server = Server(json=self.saved_json_file)
        print '''Your visualization is being rendered at
                 http://localhost:%s/
                 Visit the url in your webgl compatible browser
                 to see the animation in full glory''' % (server.port)
        server.run()

    def _display_from_ipython(self):
        # This is a hack using IPython BackgroundJobs
        # module. Once we have an IPython2.0 release
        # It can be modified to display visualizations
        # in IPython output cell itself.
        server = Server(json=self.saved_json_file)
        jobs = bg.BackgroundJobManager()
        jobs.new('server.run()')

        print '''
        Your visualization is being rendered at
        http://localhost:%s/
        Opening the visualization in new tab...'''%(server.port)
        webbrowser.open("http://localhost:%s/"%server.port)

    def display(self):
        """
        display method can be used in two ways.
        When called from IPython notebook, it shows the visualization
        in the form of output cell in the IPython notebook.
        If it is called from python interpreter or
        IPython interpreter(not notebook), It generates an html file,
        in the current directory, which can be opened in the webgl
        compliant browser for viewing the visualizations.

        The simulation data is used from this scene, hence
        all simulation data generation methods should be called before
        calling this method

        """
        try:
            # If it detects any IPython frontend
            # (qtconsole, interpreter or notebook)
            config = get_ipython().config
            self._display_from_ipython()

        except:
            self._display_from_interpreter()
