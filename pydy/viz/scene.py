#!/usr/bin/env python

# standard library
from __future__ import division
import os
import json
import distutils
import distutils.dir_util
import datetime
from collections import OrderedDict
from pkg_resources import parse_version
# external
from sympy.physics.mechanics import ReferenceFrame, Point

# local
from .camera import PerspectiveCamera
from .server import run_server
from .light import PointLight
from ..system import System

__all__ = ['Scene']

try:
    import IPython
    from IPython.html import widgets
    from IPython.display import display, Javascript

    ipy_ver = parse_version(IPython.__version__)
    if ipy_ver < parse_version('2.0'):
        raise ImportWarning('PyDy supports IPython >= 2.0.0, while you have ' +
                            'IPython ' + IPython.__version__ + ' installed. ' +
                            'IPython related functionalities will not be ' +
                            'available')
        ipython_less_than_3 = None
    elif ipy_ver >= parse_version('2.0') and ipy_ver < parse_version('3.0'):
        ipython_less_than_3 = True
    else:  # ipython >= 3.0
        ipython_less_than_3 = False

except ImportError:
    IPython = None
    ipython_less_than_3 = None


class Scene(object):
    """The Scene class holds all the data required for the visualizations
    animation of a system.

    It has methods for inputting the numerical data from the numerical
    integrations of Equations of Motions and convert them to JSON values,
    which can be then parsed by Javascripts(webgls).

    A scene object takes a ReferenceFrame, and a Point as required
    arguments. The reference_frame and point act as the inertial frame and
    origin with respect to which all objects are oriented and rendered in
    the visualizations

    A scene needs to be supplied with visualization_frames, Cameras, and
    Light objects, as optional arguments.

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
        camera : Camera, optional
            The camera with which to display the object. Default is
            PerspectiveCamera, with reference_frame and origin same as
            defined for this scene.
        """

        self.reference_frame = reference_frame
        self.origin = origin
        self.visualization_frames = list(visualization_frames)

        try:
            self.name = kwargs['name']
        except KeyError:
            self.name = 'unnamed'

        try:
            self.cameras = kwargs['cameras']
        except KeyError:
            self.cameras = [PerspectiveCamera('DefaultCamera',
                            self.reference_frame,
                            self.origin.locatenew(
                                'p_camera',
                                10*self._reference_frame.z))]

        try:
            self.lights = kwargs['lights']
        except KeyError:
            self.lights = [PointLight('DefaultLight',
                           self.reference_frame,
                           self.origin.locatenew(
                               'p_light',
                               10*self._reference_frame.z))]

    @property
    def name(self):
        """Returns the name of the scene."""
        return self._name

    @name.setter
    def name(self, new_name):
        """Sets the name of the scene."""
        if not isinstance(new_name, str):
            raise TypeError("'name' should be a valid string.")
        else:
            self._name = new_name

    @property
    def origin(self):
        """Returns the origin point of the scene."""
        return self._origin

    @origin.setter
    def origin(self, new_origin):
        """Sets the origin point of the scene."""
        if not isinstance(new_origin, Point):
            raise TypeError("'origin' should be a valid Point object.")
        else:
            self._origin = new_origin

    @property
    def reference_frame(self):
        """Returns the base reference frame of the scene."""
        return self._reference_frame

    @reference_frame.setter
    def reference_frame(self, new_reference_frame):
        """Sets the base reference frame for the scene."""
        if not isinstance(new_reference_frame, ReferenceFrame):
            raise TypeError("'reference_frame' should be a valid "
                            "ReferenceFrame object.")
        else:
            self._reference_frame = new_reference_frame

    def generate_visualization_json(self, dynamic_variables,
                                    constant_variables, dynamic_values,
                                    constant_values, fps=30,
                                    outfile_prefix=None):
        """Creates two JSON files in the current working directory. One
        contains the scene information and one contains the simulation data.

        Parameters
        ==========
        dynamic_variables : sequence of SymPy functions of time, len(m)
            The variables representing the state of the system. They should
            be in the same order as ``dynamic_values``.
        constant_variables : sequence of SymPy symbols, len(p)
            The variables representing the constants in the system. They
            should be in the same order as ``constant_variables``.
        dynamic_values : ndarray, shape(n, m)
            The trajectories of the states.
        constant_values : ndarray, shape(p,)
            The numerical values of the constants.
        fps : int, optional, default=30
            Frames per second at which animation should be displayed. Please
            not that this should not exceed the hardware limit of the
            display device to be used. Default is 30fps.
        outfile_prefix : str, optional, default=None
            A prefix for the JSON files. The files will be named as
            `outfile_prefix_scene_desc.json` and
            `outfile_prefix_simulation_data.json`. If not specified a
            timestamp shall be used as the prefix.


        """

        # TODO : The colons need to be removed from this file name.
        if outfile_prefix is None:
            timestamp = str(datetime.datetime.now())
            outfile_prefix = "_".join(timestamp.split(".")[0].split(" "))

        self.scene_json_file = outfile_prefix + "_scene_desc.json"
        self.simulation_json_file = outfile_prefix + "_simulation_data.json"

        self.constant_variables = constant_variables
        self.dynamic_variables = dynamic_variables
        self.constant_values = constant_values
        self.dynamic_values = dynamic_values
        self.outfile_prefix = outfile_prefix
        self.fps = fps

        constant_map = dict(zip(constant_variables, constant_values))
        # TODO : This assumes that all constants have unique strings and
        # that they are valid strings for the JSON file.
        constant_variables_str = [str(i) for i in constant_variables]
        constant_map_for_json = dict(zip(constant_variables_str,
                                         constant_values))

        self._simulation_data_dict = \
            self.generate_simulation_dict(dynamic_variables,
                                          constant_variables,
                                          dynamic_values,
                                          constant_values)

        self._scene_data_dict = \
            self.generate_scene_dict(constant_map=constant_map)

        self._scene_data_dict["simulationData"] = self.simulation_json_file
        # NOTE : Python 3 division is imported at the top of the file so
        # this will be a float.
        self._scene_data_dict["timeDelta"] = 1 / fps
        self._scene_data_dict["timeSteps"] = dynamic_values.shape[0]
        self._scene_data_dict["constant_map"] = constant_map_for_json

        with open(self.scene_json_file, 'w') as scene_data_outfile:
            scene_data_outfile.write(json.dumps(self._scene_data_dict,
                                                indent=4,
                                                separators=(',', ': ')))

        with open(self.simulation_json_file, 'w') as simulation_data_outfile:
            simulation_data_outfile.write(json.dumps(
                self._simulation_data_dict, indent=4,
                separators=(',', ': ')))

    def generate_scene_dict(self, constant_map={}):
        """Generates a dictionary containing all of the information needed
        to build the scene.

        Parameters
        ==========
        constant_map : dictionary
            A map of symbolic constants to numerical values. This is
            typically used if there are symbolics in the parameters of the
            Shape objects.

        Returns
        =======
        scene_info : dictionary

        """

        self._scene_info = {}
        self._scene_info["source"] = "PyDy"
        self._scene_info["name"] = self.name
        self._scene_info["newtonian_frame"] = str(self.reference_frame)
        # TODO : This should be accomodated in scene instead of width/height
        # of scene.
        self._scene_info["workspaceSize"] = 0.2

        self._scene_info["objects"] = {}
        self._scene_info["cameras"] = {}
        self._scene_info["lights"] = {}

        for frame in self.visualization_frames:
            _object_info = frame.generate_scene_dict(constant_map=constant_map)
            self._scene_info["objects"].update(_object_info)

        for camera in self.cameras:
            _object_info = camera.generate_scene_dict()
            self._scene_info["cameras"].update(_object_info)

        for light in self.lights:
            _object_info = light.generate_scene_dict()
            self._scene_info["lights"].update(_object_info)

        return self._scene_info

    def generate_simulation_dict(self, dynamic_variables,
                                 constant_variables, dynamic_values,
                                 constant_values):
        """Returns a dictionary containing all of the simulation
        information. It consists of all the simulation data along with
        references to the objects, for allowing motion to the objects in the
        PyDy visualizer.

        Parameters
        ==========
        dynamic_variables : sequence of SymPy functions of time, len(m)
            The variables representing the state of the system. They should
            be in the same order as ``dynamic_values``.
        constant_variables : sequence of SymPy symbols, len(p)
            The variables representing the constants in the system. They
            should be in the same order as ``constant_variables``.
        dynamic_values : ndarray, shape(n, m)
            The trajectories of the states.
        constant_values : ndarray, shape(p,)
            The numerical values of the constants.

        Returns
        =======
        simulation_info : dictionary

        Notes
        =====

        This method must be called before ``generate_scene_dict``.

        """
        self._simulation_info = {}

        for frame in self.visualization_frames:
            frame.generate_transformation_matrix(self.reference_frame,
                                                 self.origin)
            frame.generate_numeric_transform_function(dynamic_variables,
                                                      constant_variables)
            frame.evaluate_transformation_matrix(dynamic_values,
                                                 constant_values)

            self._simulation_info.update(frame.generate_simulation_dict())

        for camera in self.cameras:
            camera.generate_transformation_matrix(self.reference_frame,
                                                  self.origin)
            camera.generate_numeric_transform_function(dynamic_variables,
                                                       constant_variables)
            camera.evaluate_transformation_matrix(dynamic_values,
                                                  constant_values)

            self._simulation_info.update(camera.generate_simulation_dict())

        for light in self.lights:
            light.generate_transformation_matrix(self.reference_frame,
                                                 self.origin)
            light.generate_numeric_transform_function(dynamic_variables,
                                                      constant_variables)
            light.evaluate_transformation_matrix(dynamic_values,
                                                 constant_values)

            self._simulation_info.update(light.generate_simulation_dict())

        # TODO : This is bad practice. The method should either return the
        # data or mutate the object, but not both.
        return self._simulation_info

    def generate_visualization_json_system(self, system, **kwargs):
        """Creates the visualization JSON files for the provided system.

        Parameters
        ==========
        system : pydy.system.System
            A fully developed PyDy system that is prepared for the
            ``.integrate()`` method.

        Notes
        =====

        The optional keyword arguments are same as the
        ``generate_visualization_json`` method.

        """
        if not isinstance(system, System):
            # NOTE : Not sure why the next line is here. Leaving it here in
            # case it breaks something that is not tested.
            self.system = None
            msg = "{} should be a valid pydy.System object".format(system)
            raise TypeError(msg)
        else:
            self.system = system

        self.generate_visualization_json(system.states,
                                         system.constants.keys(),
                                         system.integrate(),
                                         system.constants.values(), **kwargs)

    def create_static_html(self, overwrite=False, silent=False):

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
        overwrite : boolean, optional, default=False
            If True, the directory named ``static`` in the current working
            directory will be overwritten.
        Silent : boolean, optional, default=False
            If True, no messages will be displayed to
            STDOUT
        """

        dst = os.path.join(os.getcwd(), 'static')

        if os.path.exists(dst) and overwrite is False:
            ans = raw_input("The 'static' directory already exists. Would "
                            + "you like to overwrite the contents? [y|n]\n")
            if ans == 'y':
                distutils.dir_util.remove_tree(dst)
            else:
                if not silent: print "Aborted!"
                return

        src = os.path.join(os.path.dirname(__file__), 'static')
        if not silent: print("Copying static data.")
        distutils.dir_util.copy_tree(src, dst)
        if not silent: print("Copying Simulation data.")
        _scene_outfile_loc = os.path.join(os.getcwd(), 'static', self.scene_json_file)
        _simulation_outfile_loc = os.path.join(os.getcwd(), 'static', self.simulation_json_file)
        scene_outfile = open(_scene_outfile_loc, "w")
        simulation_outfile = open(_simulation_outfile_loc, "w")

        scene_outfile.write(json.dumps(self._scene_data_dict, indent=4,
                                separators=(',', ': ')))
        scene_outfile.close()
        simulation_outfile.write(json.dumps(self._simulation_data_dict, indent=4,
                                separators=(',', ': ')))
        simulation_outfile.close()
        if not silent:
            print("To view the visualization, open {}".format(
                  os.path.join(dst, 'index.html')) +
                  " in a WebGL compliant browser.")

    def remove_static_html(self, force=False):
        """Removes the ``static`` directory from the current working
        directory.

        Parameters
        ----------
        force : boolean, optional, default=False
            If true, no warning is issued before the removal of the
            directory.

        """
        if not os.path.exists('static'):
            print "All Done!"
            return

        if not force:
            ans = raw_input("Are you sure you would like to delete the " +
                            "'static' directory? [y|n]\n")
            if ans == 'y':
                force = True

        if force:
            distutils.dir_util.remove_tree(os.path.join(os.getcwd(),
                                                        'static'))
            print "All Done!"
        else:
            print "aborted!"

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
        self.create_static_html()
        run_server(scene_file=self.scene_json_file)

    def display_ipython(self):
        """
        Method to display the visualization inside the
        IPython notebook. It is only supported by IPython
        versions>=2.0.0

        """

        # Raise error whenever display_ipython() is called
        # and IPython is not installed or IPython < '2.0.0'
        if IPython is None:
            raise ImportError('IPython is not installed but is required. ' +
                              'Please installed IPython >= 2.0 and try again')
        elif ipython_less_than_3 is None:
            raise ImportError('You have IPython ' + IPython.__version__ +
                              ' installed but PyDy supports IPython >= 2.0.0' +
                              'Please update IPython and try again')

        self.create_static_html(silent=True)
        self._widget_dict = OrderedDict()
        if ipython_less_than_3:
            self.container = widgets.ContainerWidget()
        else:
            self.container = widgets.Box()
        components = []
        for var, init_val in \
            zip(self.constant_variables, self.constant_values):
            if ipython_less_than_3:
                self._widget_dict[str(var)] = widgets.FloatTextWidget(
                    value=init_val,
                    description=str(var))
            else:
                self._widget_dict[str(var)] = widgets.FloatText(
                    value=init_val,
                    description=str(var))
            components.append(self._widget_dict[str(var)])

        if ipython_less_than_3:
            self.button = widgets.ButtonWidget(description="Rerun Simulations")
        else:
            self.button = widgets.Button(description="Rerun Simulations")

        def button_click(clicked):

            if ipython_less_than_3:
                self.button.add_class('disabled')
            else:
                self.button._dom_classes = ['disabled']

            self.button.description = 'Rerunning Simulation ...'
            self.constant_values = []

            for i in self._widget_dict.values():
                self.constant_values.append(i.value)

            if self.system is not None:
                # update system constants
                self.system.constants = dict(zip(self.system.constants,
                                                 self.constant_values))
                self.generate_visualization_json_system(self.system)
            else:
                self.generate_visualization_json(
                    self.dynamic_variables,
                    self.constant_variables, self.dynamic_values,
                    self.constant_values,
                    fps=self.fps,
                    outfile_prefix=self.outfile_prefix)

            self.create_static_html(overwrite=True, silent=True)
            js = 'jQuery("#json-input").val("{}");'.format('static/' +
                                                           self.scene_json_file)
            display(Javascript(js))
            display(Javascript('jQuery("#simulation-load").click()'))

            if ipython_less_than_3:
                self.button.remove_class('disabled')
            else:
                self.button._dom_classes = ['enabled']

            self.button.description = 'Rerun Simulation'

        self.button.on_click(button_click)
        html_file = open("static/index_ipython.html")

        if ipython_less_than_3:
            self.html_widget = widgets.HTMLWidget(
                value=html_file.read().format(load_url='static/' +
                                              self.scene_json_file))
        else:
            self.html_widget = widgets.HTML(
                value=html_file.read().format(load_url='static/' +
                                              self.scene_json_file))
        self.container.children = components

        if ipython_less_than_3:
            self.container.set_css({"max-height": "10em",
                                    "overflow-y": "scroll",
                                    "display": "block"
                                    })
            self.html_widget.set_css({"display": "block",
                                      "float": "left"
                                      })
        else:
            self.container._css = [("canvas", "width", "100%")]

        display(self.container)
        display(self.button)
        display(self.html_widget)

        if ipython_less_than_3:
            self.button.add_class('btn-info')
        else:
            self.button._dom_classes = ['btn-info']
