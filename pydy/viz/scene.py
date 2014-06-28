#!/usr/bin/env python

# standard library
import os
import json
import shutil
import distutils
import distutils.dir_util
import webbrowser

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

    def generate_visualization_dict(self, dynamic_variables,
                                    constant_variables, dynamic_values,
                                    constant_values):
        """
        generate_visualization_dict() method generates
        a dictionary of visualization data


        Parameters
        ==========
        dynamic_variables : sequence of sympy.Functions
            This sequence contains all the functions of time which are
            required for generating the transformation matrices of all the
            visualization frames in the scene.
        constant_variables : Sympifyable list or tuple
            This contains all the symbols for the parameters which are
            used for defining various objects in the system.
        dynamic_values : list or tuple
            initial states of the system. The list or tuple
            should be respective to the state_sym.
        constant_values : list or tuple
            values of the parameters. The list or tuple
            should be respective to the par_sym.

        Returns
        =======

        The dictionary contains following keys:
        1) Width of the scene.
        2) Height of the scene.
        3) name of the scene.
        4) frames in the scene, which contains sub-dictionaries
           of all the visualization frames information.


        """

        self._scene_data = {}
        self._scene_data['name'] = self._name
        self._scene_data['height'] = self._height
        self._scene_data['width'] = self._width
        self._scene_data['frames'] = []
        self._scene_data['cameras'] = []
        self._scene_data['lights'] = []

        constant_map = dict(zip(constant_variables, constant_values))

        for frame in self.visualization_frames:
            frame.generate_transformation_matrix(self._reference_frame,
                                                 self._origin)
            frame.generate_numeric_transform_function(dynamic_variables,
                                                      constant_variables)
            frame.evaluate_transformation_matrix(
                dynamic_values, constant_values)

            self._scene_data['frames'].append(
                frame.generate_visualization_dict(constant_map=constant_map))

        for camera in self.cameras:
            camera.generate_transformation_matrix(self._reference_frame,
                                                  self._origin)
            camera.generate_numeric_transform_function(dynamic_variables,
                                                       constant_variables)
            camera.evaluate_transformation_matrix(dynamic_values,
                                                  constant_values)

            self._scene_data['cameras'].append(
                camera.generate_visualization_dict()
            )

        for light in self.lights:
            light.generate_transformation_matrix(self._reference_frame,
                                                 self._origin)
            light.generate_numeric_transform_function(dynamic_variables,
                                                      constant_variables)
            light.evaluate_transformation_matrix(dynamic_values,
                                                 constant_values)

            self._scene_data['lights'].append(
                light.generate_visualization_dict())

        return self._scene_data


    def generate_visualization_json(self, dynamic_variables,
                                    constant_variables, dynamic_values,
                                    constant_values, save_to='data.json'):
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

        save_to : str
            path to the file where to write the generated data JSON.
            the path should be chosen such as to have the write
            permissions to the user.

        Returns
        =======

        The dictionary contains following keys:
        1) Width of the scene.
        2) Height of the scene.
        3) name of the scene.
        4) frames in the scene, which contains sub-dictionaries
           of all the visualization frames information.


        """
        self.saved_json_file = save_to
        self._data_dict = self.generate_visualization_dict(dynamic_variables,
                                                           constant_variables,
                                                           dynamic_values,
                                                           constant_values)

        outfile = open(self.saved_json_file, 'w')
        outfile.write(json.dumps(self._data_dict, indent=4,
                                 separators=(',', ': ')))
        outfile.close()

    def create_static_html(self, overwrite=False):
        """Creates a directory named ``static`` in the current working
        directory which contains all of the HTML, CSS, and Javascript files
        required to run the visualization. Simply open ``static/index.html``
        in a WebGL compliant browser to view and interact with the
        visualization.

        This method can also be used to output files for embedding the
        visualizations in the static webpages. Simply copy the contents of
        static directory in the relevant directory for embedding in a static
        website.

        Parameters
        ----------
        dir_name : string
            A valid directory name.
        overwrite : boolean, optional, default=False
            If true, the directory named ``static`` in the current working
            directory will be overwritten.

        """

        dst = os.path.join(os.getcwd(), 'static')

        if os.path.exists(dst) and overwrite is False:
            ans = raw_input("The 'static' directory already exists. Would "
                            + "you like to overwrite the contents? [y|n]\n")
            if ans == 'y':
                shutil.rmtree(dst)
                overwrite = True
        elif os.path.exists(dst) and overwrite is True:
            shutil.rmtree(dst)
        elif not os.path.exists(dst):
            overwrite = True

        if overwrite is True:
            src = os.path.join(os.path.dirname(__file__), 'static')
            print("Copying static data.")
            shutil.copytree(src, dst)
            print("Copying Simulation data.")
            _outfile_loc = os.path.join(os.getcwd(), 'static', 'data.json')
            outfile = open(_outfile_loc, "w")
            # For static rendering, we need to define json data as a
            # JavaScript variable.
            outfile.write('var JSONObj=')
            outfile.write(json.dumps(self._data_dict, indent=4,
                                    separators=(',', ': ')))
            outfile.write(';')
            outfile.close()
            print("To view the visualization, open {}".format(
                os.path.join(dst, 'index.html')) +
                " in a WebGL compliant browser.")
        else:
            print('Aborted.')

    def remove_static_html(self, force=False):
        """Removes the ``static`` directory from the current working
        directory.

        Parameters
        ----------
        force : boolean, optional, default=False
            If true, no warning is issued before the removal of the
            directory.

        """
        if os.path.exists('static'):
            if force is False:
                ans = raw_input("Are you sure you would like to delete the " +
                                "'static' directory? [y|n]\n")
                if ans == 'y':
                    force = True

            if force is True:
                print 'Cleaning up static directory..'
                distutils.dir_util.remove_tree(os.path.join(os.getcwd(),
                                                            'static'))
                print 'All Done!'
            else:
                print('Aborted.')

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
