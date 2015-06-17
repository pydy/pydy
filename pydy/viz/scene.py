#!/usr/bin/env python

# standard library
from __future__ import division
import os
import sys
import warnings
import json
import distutils
import distutils.dir_util
import datetime
from collections import OrderedDict
from pkg_resources import parse_version

# external
import numpy as np
from sympy import latex
from sympy.physics.mechanics import ReferenceFrame, Point

# local
from .camera import PerspectiveCamera
from .server import run_server
from .light import PointLight
from ..system import System
from ..utils import PyDyImportWarning, PyDyDeprecationWarning

if sys.version_info > (3, 0):
    raw_input = input

__all__ = ['Scene']

warnings.simplefilter('once', PyDyImportWarning)
warnings.simplefilter('once', PyDyDeprecationWarning)

try:
    import IPython
    if parse_version(IPython.__version__) < parse_version('3.0'):
        msg = ('PyDy only supports IPython >= 3.0.0. You have '
               'IPython {} installed. IPython related functionalities will '
               'not be available')
        warnings.warn(msg.format(IPython.__version__), PyDyImportWarning)
        ipython_less_than_3 = True
    else:
        ipython_less_than_3 = False
        from IPython.html import widgets
        from IPython.display import display, Javascript
except ImportError:
    IPython = None


class Scene(object):
    """The Scene class generates all of the data required for animating a
    set of visualization frames.

    """
    def __init__(self, reference_frame, origin, *visualization_frames,
                 **kwargs):
        """Initialize a Scene instance.

        Parameters
        ==========
        reference_frame : sympy.physics.mechanics.ReferenceFrame
            The base reference frame for the scene. The motion of all of the
            visualization frames, cameras, and lights will be generated with
            respect to this reference frame.
        origin : sympy.physics.mechanics.Point
            The base point for the scene. The motion of all of the
            visualization frames, cameras, and lights will be generated with
            respect to this point.
        visualization_frames : VisualizationFrame
            One or more visualization frames which are to be displayed in
            the scene.
        name : string, optional, default='unnamed'
            Name of Scene object.
        cameras : list of Camera instances, optional
            The cameras with which to display the object. The first camera
            is used to display the scene initially. The default is a single
            PerspectiveCamera tied to the base reference frame and
            positioned away from the origin along the reference frame's z
            axis.
        lights : list of Light instances, optional
            The lights used in the scene. The default is a single Light tied
            to the base reference frame and positioned away from the origin
            along the reference frame's z axis at the same point as the
            default camera.
        system : System, optional, default=None
            A PyDy system class which is initiated such that the
            ``integrate()`` method will produce valid state trajectories.
        time : array_like, shape(n,), optional, default=None
            Monotoncially increaing float values of time that correspond to
            the state trajectories.
        constants : dictionary, optional, default=None
            A dictionary that maps SymPy symbols to floats. This should
            contain at least all necessary symbols to evaluate the
            transformation matrices of the visualization frame, cameras, and
            lights and to evaluate the Shapes' parameters.
        states_symbols : sequence of functions of time, len(m), optional, default=None
            An ordered sequence of the SymPy functions that represent the
            states. The order must match the order of the
            ``states_trajectories``.
        states_trajectories : array_like, shape(n, m), optional, default=None
            A two dimensional array with numerical values for each state at
            each poitn in time during the animation.

        Notes
        =====
        The user is allow to supply either system or time, constants,
        states_symbols, and states_trajectories. Providing a System allows for
        interactively changing the simulation parameters via the Scene GUI
        in the IPython notebook.

        """

        self.reference_frame = reference_frame
        self.origin = origin
        self.visualization_frames = list(visualization_frames)

        vec = 10 * self.reference_frame.z
        self._default_camera_point = self.origin.locatenew('p_camera', vec)
        self._default_light_point = self.origin.locatenew('p_light', vec)

        default_kwargs = {'name': 'unnamed',
                          'cameras': [PerspectiveCamera('DefaultCamera',
                                                        self.reference_frame,
                                                        self._default_camera_point)],
                          'lights': [PointLight('DefaultLight', self.reference_frame,
                                                self._default_light_point)],
                          'system': None,
                          'time': None,
                          'constants': None,
                          'states_symbols': None,
                          'states_trajectories': None,
                          'frames_per_second': 30}
        default_kwargs.update(kwargs)
        for k, v in default_kwargs.items():
            setattr(self, k, v)

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

    @property
    def system(self):
        return self._system

    @system.setter
    def system(self, new_system):

        if new_system is not None and not isinstance(new_system, System):
            msg = "{} should be a valid pydy.System object".format(new_system)
            raise TypeError(msg)

        if new_system is not None:
            msg = ('The {} attribute has already been set, so the system '
                   'cannot be set. Use the clear_trajectories method to '
                   'set all relevant attributes to None.')
            for attr in ['time', 'constants', 'states_symbols',
                         'states_trajectories']:
                try:
                    if getattr(self, attr) is not None:
                        raise ValueError(msg.format(attr))
                except AttributeError:
                    pass

        self._system = new_system

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, new_time):

        try:
            if new_time is not None and self.system is not None:
                msg = ('The system attribute has already been set, so the '
                       'time cannot be set. Set Scene.system = None to '
                       'allow a time array to be added.')
                raise ValueError(msg)
        except AttributeError:
            pass

        try:
            if new_time is not None and self.states_trajectories is not None:
                len_traj = self.states_trajectories.shape[0]
                if len(new_time) != len_traj:
                    msg = ('The time array length, {}, does not match the '
                           'length of the state trajectories array, {}.')
                    raise ValueError(msg.format(len(new_time), len_traj))
        except AttributeError:
            pass

        if new_time is None:
            self._time = new_time
        else:
            self._time = np.asarray(new_time)

    @property
    def states_symbols(self):
        return self._states_symbols

    @states_symbols.setter
    def states_symbols(self, new_states_symbols):

        try:
            if new_states_symbols is not None and self.system is not None:
                msg = ('The system attribute has already been set, so the '
                       'coordinates cannot be set. Set Scene.system = None '
                       'to allow a coordinates array to be added.')
                raise ValueError(msg)
        except AttributeError:
            pass

        try:
            if (new_states_symbols is not None and self.states_trajectories
                is not None):
                len_traj = self.states_trajectories.shape[1]
                if len(new_states_symbols) != len_traj:
                    msg = ('The number of states, {}, does not match the '
                           'number of states present in the state '
                           'trajectories array, {}.')
                    raise ValueError(msg.format(len(new_states_symbols), len_traj))
        except AttributeError:
            pass

        self._states_symbols = new_states_symbols

    @property
    def states_trajectories(self):
        return self._states_trajectories

    @states_trajectories.setter
    def states_trajectories(self, new_states_trajectories):

        try:
            if new_states_trajectories is not None and self.system is not None:
                msg = ('The system attribute has already been set, so the '
                       'states_trajectories cannot be set. Set Scene.system '
                       '= None to allow a states_trajectories array to be '
                       'added.')
                raise ValueError(msg)
        except AttributeError:
            pass

        try:
            if new_states_trajectories is not None and self.time is not None:
                if len(self.time) != new_states_trajectories.shape[0]:
                    msg = ("The number of time instances do not match the "
                           "number in the time array.")
                    raise ValueError(msg)
        except AttributeError:
            pass

        try:
            if (new_states_trajectories is not None and self.states_symbols
                is not None):
                if new_states_trajectories.shape[1] != len(self.states_symbols):
                    msg = ("The number of states in the trajectory do not "
                           "match the number of states symbols.")
                    raise ValueError(msg)
        except AttributeError:
            pass

        self._states_trajectories = new_states_trajectories

    @property
    def constants(self):
        return self._constants

    @constants.setter
    def constants(self, new_constants):
        try:
            if new_constants is not None and self.system is not None:
                msg = ('The system attribute has already been set, so the '
                       'constants cannot be set. Set Scene.system = None to '
                       'allow a constants array to be added.')
                raise ValueError(msg)
        except AttributeError:
            pass
        self._constants = new_constants

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
        warnings.warn("This method will be removed in PyDy 0.4.0, set these "
                      "values through the proper attributes instead.",
                      PyDyDeprecationWarning)

        self.states_symbols = dynamic_variables
        self.states_trajectories = dynamic_values
        self.constants = dict(zip(constant_variables, constant_values))
        self.frames_per_second = fps
        self._generate_json(prefix=outfile_prefix)

    def _generate_json(self, directory=None, prefix=None):
        """Creates two JSON files in the specified directory. One contains
        the scene information and one contains the simulation data. If
        ``directory`` is None the files will be created in the current
        working directory."""

        if directory is None:
            directory = os.getcwd()

        if prefix is None:
            prefix = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

        self._scene_json_file = prefix + "_scene_desc.json"
        self._simulation_json_file = prefix + "_simulation_data.json"

        if self.system is None:
            constants = self.constants
            num_time_steps = self.states_trajectories.shape[0]
        else:
            constants = self.system.constants
            num_time_steps = len(self.system.times)
        # TODO : This assumes that all constants have unique strings and
        # that they are valid strings for the JSON file which can possibly
        # fail.
        constant_map_for_json = {str(k): v for k, v in constants.items()}

        self._generate_simulation_dict()
        self._generate_scene_dict()

        self._scene_info["simulationData"] = self._simulation_json_file
        # NOTE : Python 3 division is imported at the top of the file so
        # this will be a float.
        self._scene_info["timeDelta"] = 1 / self.frames_per_second
        self._scene_info["timeSteps"] = num_time_steps
        self._scene_info["constant_map"] = constant_map_for_json

        scene_file_path = os.path.join(directory, self._scene_json_file)
        with open(scene_file_path, 'w') as scene_data_outfile:
            scene_data_outfile.write(json.dumps(self._scene_info,
                                                indent=4,
                                                separators=(',', ': ')))

        sim_file_path = os.path.join(directory, self._simulation_json_file)
        with open(sim_file_path, 'w') as simulation_data_outfile:
            simulation_data_outfile.write(json.dumps(
                self._simulation_info, indent=4,
                separators=(',', ': ')))

    def _generate_scene_dict(self):
        """Generates a dictionary containing all of the information needed
        to build the scene."""

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
            if self.system is None:
                constants = self.constants
            else:
                constants = self.system.constants
            object_info = frame.generate_scene_dict(constant_map=constants)
            self._scene_info["objects"].update(object_info)

        for camera in self.cameras:
            object_info = camera.generate_scene_dict()
            self._scene_info["cameras"].update(object_info)

        for light in self.lights:
            object_info = light.generate_scene_dict()
            self._scene_info["lights"].update(object_info)

    def _generate_simulation_dict(self):
        """Returns a dictionary containing all of the simulation
        information. It consists of all the simulation data along with
        references to the objects, for allowing motion to the objects in the
        PyDy visualizer.

        Notes
        =====

        This method must be called before ``generate_scene_dict``.

        """
        if self.system is None:
            constants_symbols = self.constants.keys()
            constants_values = self.constants.values()
            states_symbols = self.states_symbols
            states_trajectories = self.states_trajectories
        else:
            constants_symbols = self.system.constants.keys()
            constants_values = self.system.constants.values()
            states_symbols = self.system.states
            states_trajectories = self.system.integrate()

        self._simulation_info = {}

        for group in [self.visualization_frames, self.cameras, self.lights]:
            for frame in group:
                frame.generate_transformation_matrix(self.reference_frame,
                                                     self.origin)
                frame.generate_numeric_transform_function(states_symbols,
                                                          constants_symbols)
                frame.evaluate_transformation_matrix(states_trajectories,
                                                     constants_values)
                self._simulation_info.update(frame.generate_simulation_dict())

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

        self.system = system
        prefix = None
        for k, v in kwargs.items():
            if k == 'fps':
                self.frames_per_second = v
            if k == 'outfile_prefix':
                prefix = v
        self._generate_json(prefix=prefix)

    def create_static_html(self, overwrite=False, silent=False):
        """Creates a directory named ``static`` in the current working
        directory which contains all of the HTML, CSS, and Javascript files
        required to run the vizualization application. To run the
        application, navigate into the ``static`` directory and start a
        webserver from that directory, e.g.::

            $ python -m SimpleHTTPServer

        Now, in a WebGL compliant browser, navigate to::

            http://127.0.0.1:8000

        to view and interact with the visualization.

        This method can also be used to output files for embedding the
        visualizations in static webpages. Simply copy the contents of
        static directory in the relevant directory for embedding in a static
        website.

        Parameters
        ----------
        overwrite : boolean, optional, default=False
            If True, the directory named ``static`` in the current working
            directory will be overwritten.
        silent : boolean, optional, default=False
            If True, no messages will be displayed to STDOUT.

        """

        dst = os.path.join(os.getcwd(), 'static')

        if os.path.exists(dst) and overwrite is False:
            ans = raw_input("The 'static' directory already exists. Would "
                            "you like to overwrite the contents? [y|n]\n")
            if ans == 'y':
                distutils.dir_util.remove_tree(dst)
            else:
                if not silent:
                    print("Aborted!")
                return

        src = os.path.join(os.path.dirname(__file__), 'static')

        if not silent:
            print("Copying static data.")
        distutils.dir_util.copy_tree(src, dst)

        if not silent:
            print("Copying Simulation data.")

        _scene_outfile_loc = os.path.join(os.getcwd(), 'static',
                                          self._scene_json_file)
        _simulation_outfile_loc = os.path.join(os.getcwd(), 'static',
                                               self._simulation_json_file)
        scene_outfile = open(_scene_outfile_loc, "w")
        simulation_outfile = open(_simulation_outfile_loc, "w")

        scene_outfile.write(json.dumps(self._scene_info, indent=4,
                                       separators=(',', ': ')))
        scene_outfile.close()
        simulation_outfile.write(json.dumps(self._simulation_info,
                                            indent=4,
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
            print("All Done!")
            return

        if not force:
            ans = raw_input("Are you sure you would like to delete the " +
                            "'static' directory? [y|n]\n")
            if ans == 'y':
                force = True

        if force:
            distutils.dir_util.remove_tree(os.path.join(os.getcwd(),
                                                        'static'))
            print("All Done!")
        else:
            print("aborted!")

    def display(self):
        """Displays the scene in the default webbrowser."""
        self.create_static_html()
        run_server(scene_file=self._scene_json_file)

    def _rerun_button_callback(self, btn):
        """Callback for the "Rerun Simulation" button. When executed the
        parameter values are collected from the text input widgets and used
        in a new simulation of the model."""

        btn._dom_classes = ['btn-info', 'active', 'disabled']

        btn.description = 'Rerunning Simulation...'

        original_scene_file = self._scene_json_file
        original_constants = self._system.constants
        try:
            self._system.constants = {s: w.value for s, w in
                                      self._constants_text_widgets.items()}
            self.generate_visualization_json_system(self._system)
        except:
            print('Simulation rerun failed, using previous simulation data.')
            # If the simulation fails for any reason we revert everything
            # back to the previous state, including filling the text widgets
            # with the previous values of the constants. The _scene_info and
            # _simulation_data dicts may be in a bad state, but that should
            # be ok, because generate_visualiation_json will have to be run
            # again for anything new to happen.
            self._scene_json_file = original_scene_file
            self._system.constants = original_constants
            self._fill_constants_widgets()

        self.create_static_html(overwrite=True, silent=True)

        js_tmp = 'jQuery("#json-input").val("{}");'
        js = js_tmp.format('static/' + self._scene_json_file)
        display(Javascript(js))
        display(Javascript('jQuery("#simulation-load").click()'))

        btn._dom_classes = ['btn-info', 'enabled']

        btn.description = self._rerun_button_desc

    def _fill_constants_widgets(self):
        """Fills up the constants widget with the current constants symbols
        and values."""

        for sym, init_val in self._system.constants.items():

            desc = latex(sym, mode='inline')

            text_widget = widgets.FloatText(value=init_val,
                                            description=desc)

            self._constants_text_widgets[sym] = text_widget

    def _initialize_rerun_button(self):
        """Construct a button for controlling rerunning the simulations."""

        self._rerun_button = widgets.Button()
        self._rerun_button._dom_classes = ['btn-info']

        self._rerun_button_desc = "Rerun Simulation"
        self._rerun_button.description = self._rerun_button_desc

        self._rerun_button.on_click(self._rerun_button_callback)

    def display_ipython(self):
        """Displays the scene using an IPython widget inside an IPython
        notebook cell.

        Notes
        =====
        IPython widgets are only supported by IPython versions >= 3.0.0.

        """

        if IPython is None:
            raise ImportError('IPython is not installed but is required. ' +
                              'Please install IPython >= 3.0 and try again')

        if ipython_less_than_3:
            msg = ('You have IPython {} installed but PyDy only supports '
                   'IPython >= 3.0. Please update IPython and try again.')
            raise ImportError(msg.format(IPython.__version__))

        self.create_static_html(silent=True)

        # Only create the constants input boxes and the rerun simulation
        # button if the scene was generated with a System.
        if self._system is not None:

            # Construct a container that holds all of the constants input
            # text widgets.
            self._constants_container = widgets.Box()
            self._constants_container._css = [("canvas", "width", "100%")]

            self._constants_text_widgets = OrderedDict()
            self._fill_constants_widgets()
            # Add all of the constants widgets to the container.
            self._constants_container.children = \
                tuple(v for v in self._constants_text_widgets.values())

            self._initialize_rerun_button()

            display(self._constants_container)
            display(self._rerun_button)

        with open("static/index_ipython.html", 'r') as html_file:
            html = html_file.read()

        html = html.format(load_url='static/' + self._scene_json_file)

        self._html_widget = widgets.HTML(value=html)

        display(self._html_widget)
