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
from math import sqrt

# external
from pkg_resources import parse_version
import numpy as np
from sympy import latex
from sympy.physics.mechanics import ReferenceFrame, Point
import pythreejs as p3js

# local
from .camera import PerspectiveCamera
from .server import Server
from .light import PointLight
from .trajectory_link import trajectory_link, play_link
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
        # If IPython >= 4.0 use ipywidgets if installed to avoid the
        # deprecation warning.
        try:
            import ipywidgets as widgets
        except ImportError:
            from IPython.html import widgets
        from IPython.display import display, Javascript
        from IPython.utils.traitlets import CFloat, List
except ImportError:
    IPython = None


def quaternion_from_rotation(m):
    """
    m is a 3 by 3 matrix, as a list of rows.  The columns of this matrix are
    the vectors x, y, and z
    """
    x = m[0:3]
    y = m[3:6]
    z = m[6:9]
    trace = x[0]+y[1]+z[2]
    if (trace>0):
        s = 0.5/sqrt(trace+1)
        quaternion = [(y[2]-z[1])*s, (z[0]-x[2])*s, (x[1]-y[0])*s, 0.25/s]
    elif (x[0]>y[1] and x[0]>z[2]):
        s = 2.0*sqrt(1.0+x[0]-y[1]-z[2])
        quaternion = [0.25*s, (y[0]+x[1])/s, (z[0]+x[2])/s, (y[2]-z[1])/s]
    elif (y[1]>z[2]):
        s = 2.0*sqrt(1.0+y[1]-x[0]-z[2])
        quaternion = [(y[0]+x[1])/s, 0.25*s, (z[1]+y[2])/s, (z[0]-x[2])/s]
    else:
        s = 2.0*sqrt(1.0+z[2]-x[0]-y[1])
        quaternion = [(z[0]+x[2])/s, (z[1]+y[2])/s, 0.25*s, (x[1]-y[0])/s]
    return quaternion


def _create_pythreejs_geometery(shape):
    from . import shapes
    if type(shape) == shapes.Cube:
        geo = p3js.BoxGeometry(length=shape.length)
    elif type(shape) == shapes.Cylinder:
        geo = p3js.CylinderGeometry(radiusTop=shape.radius,
                                    radiusBottom=shape.radius,
                                    height=shape.length)
    elif type(shape) == shapes.Cone:
        geo = p3js.CylinderGeometry(radiusTop=0,
                                    radiusBottom=shape.radius,
                                    height=shape.length)
    elif type(shape) == shapes.Sphere:
        geo = p3js.SphereGeometry(radius=shape.radius)
    elif type(shape) == shapes.Circle:
        geo = p3js.CircleGeometry(radius=shape.radius, segments=32)
    elif type(shape) == shapes.Plane:
        geo = p3js.PlaneGeometry(width=shape.width, length=shape.length)
    elif type(shape) == shapes.Tetrahedron:
        geo = p3js.TetrahedronGeometry(radius=shape.radius)
    elif type(shape) == shapes.Octahedron:
        geo = p3js.OctahedronGeometry(radius=shape.radius)
    elif type(shape) == shapes.Icosahedron:
        geo = p3js.IcosahedronGeometry(radius=shape.radius)
    elif type(shape) == shapes.Torus:
        geo = p3js.TorusGeometry(radius=shape.radius, tube=shape.tube_radius)
    elif type(shape) == shapes.TorusKnot:
        geo = p3js.TorusKnotGeometry(radius=shape.radius,
                                     tube=shape.tube_radius)
    else:
        raise TypeError('No pythreejs geometry for {}'.format(type(shape)))
    return geo


class Scene(object):
    """The Scene class generates all of the data required for animating a
    set of visualization frames.

    """

    pydy_directory = "pydy-resources"

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
        times : array_like, shape(n,), optional, default=None
            Monotoncially increaing float values of time that correspond to
            the state trajectories.
        constants : dictionary, optional, default=None
            A dictionary that maps SymPy symbols to floats. This should
            contain at least all necessary symbols to evaluate the
            transformation matrices of the visualization frame, cameras, and
            lights and to evaluate the Shapes' parameters.
        states_symbols : sequence of functions, len(m), optional, default=None
            An ordered sequence of the SymPy functions that represent the
            states. The order must match the order of the
            ``states_trajectories``.
        states_trajectories : array_like, shape(n, m), optional, default=None
            A two dimensional array with numerical values for each state at
            each point in time during the animation.

        Notes
        =====
        The user is allowed to supply either system or times, constants,
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
                          'lights': [PointLight('DefaultLight',
                                                self.reference_frame,
                                                self._default_light_point)],
                          'system': None,
                          'times': None,
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
            for attr in ['times', 'constants', 'states_symbols',
                         'states_trajectories']:
                try:
                    if getattr(self, attr) is not None:
                        raise ValueError(msg.format(attr))
                except AttributeError:
                    pass

        self._system = new_system

    @property
    def times(self):
        return self._times

    @times.setter
    def times(self, new_times):

        try:
            if new_times is not None and self.system is not None:
                msg = ('The system attribute has already been set, so the '
                       'times cannot be set. Set Scene.system = None to '
                       'allow a time array to be added.')
                raise ValueError(msg)
        except AttributeError:
            pass

        try:
            if new_times is not None and self.states_trajectories is not None:
                len_traj = self.states_trajectories.shape[0]
                if len(new_times) != len_traj:
                    msg = ('The times array length, {}, does not match the '
                           'length of the state trajectories array, {}.')
                    raise ValueError(msg.format(len(new_times), len_traj))
        except AttributeError:
            pass

        if new_times is None:
            self._times = new_times
        else:
            self._times = np.array(new_times)

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
            if new_states_trajectories is not None and self.times is not None:
                if len(self.times) != new_states_trajectories.shape[0]:
                    msg = ("The number of time instances do not match the "
                           "number in the times array.")
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

    def clear_trajectories(self):
        """Sets the 'system', 'times', 'constants', 'states_symbols', and
        'states_trajectories' to None."""
        for attr in ['system', 'times', 'constants', 'states_symbols',
                     'states_trajectories']:
            setattr(self, attr, None)

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
        """Creates two JSON files and copies all the necessary static files
        in the specified directory. One of the JSON files contains the scene
        information and other one contains the simulation data.

        Parameters
        ==========
        directory : string, optional, default=None
            The directory in which the json files are placed. If None, this
            will be the current working directory.
        prefix : string
            A custom prefix for the two json files. If None, a time stamp is
            used.

        """

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

        times = None
        if self.times is not None:
            times = self.times
        elif self.system and self.system.times is not None:
            times = self.system.times

        if times is None:
            self._scene_info["timeDelta"] = 1.0 / self.frames_per_second
            self._scene_info["startTime"] = 0.0
        else:
            # Assume that times is evenly spaced and monotonic.
            # TODO: Interpolate if times are not evenly spaced.
            total_time = times[-1] - times[0]
            self._scene_info["timeDelta"] = total_time / (num_time_steps - 1)
            self._scene_info["startTime"] = times[0]

        self._scene_info["fps"] = self.frames_per_second
        self._scene_info["speedup"] = (self.frames_per_second *
                                       self._scene_info["timeDelta"])
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
                constants = self._system.constants
            object_info = frame.generate_scene_dict(constant_map=constants)
            self._scene_info["objects"].update(object_info)

        for camera in self.cameras:
            object_info = camera.generate_scene_dict()
            self._scene_info["cameras"].update(object_info)

        for light in self.lights:
            object_info = light.generate_scene_dict()
            self._scene_info["lights"].update(object_info)

    def _generate_constants_widget(self):
        """Generates an ipywidget.FlexBox containing widgets for all system
        constants, if any."""
        constants_widgets = []
        for k, v in self._system.constants.items():
            w = widgets.FloatText(value=v, width = 80,
                                  description=latex(k, mode='inline'))
            w._symbol = k
            constants_widgets.append(w)
        self._constants_widget = widgets.HBox(
                constants_widgets,
                overflow_x='scroll', width='100%')

    def _generate_time_widget(self):
        """Generates an ipywidget.FlexBox containing widgets for simulation time
        parameters."""
        times = None
        if self.times is not None:
            times = self.times
        elif self._system and self._system.times is not None:
            times = self._system.times

        self._time_widget = widgets.HBox(width='100%')
        if times is not None:
            t0 = times[0]
            tf = times[-1]
            #TODO: assumes 'times' is evenly spaced.
            dt = (tf - t0) / (len(times) - 1)
            width = 80
            self._time_widget.children = (
                widgets.FloatText(value=t0, width=width,
                                  description='start time'),
                widgets.FloatText(value=tf, width=width,
                                  description='end time'),
                widgets.FloatText(value=dt, width=width,
                                  description='time step'))

    def _generate_play_widget(self):
        play = widgets.ToggleButton(description='play')
        speedup = widgets.FloatText(description='speedup', value=1.0, width=80)
        loop = widgets.Checkbox(description='loop')
        resim = widgets.Button(description='resimulate')
        proxy = widgets.Button(description='proxy', width=5)
        controls = widgets.HBox()
        controls.children = (play, loop, proxy, speedup, resim)

        if self._time_widget.children:
            t0 = self._time_widget.children[0].value
            tf = self._time_widget.children[1].value
            dt = self._time_widget.children[2].value
            slider_type = widgets.FloatSlider
            desc = 'time:'
        else:
            t0 = 0
            dt = 1
            tf = self.states_trajectories.shape[0] - 1
            slider_type = widgets.IntSlider
            desc = 'step:'
        time_slider = slider_type(value=t0, min=t0, max=tf, step=dt,
                                  description=desc)

        self._play_widget = widgets.VBox()
        self._play_widget.children = (controls, time_slider)
        self._play_widget.border_width = 5
        self._play_widget.border_color = 'white'
        self._play_link = play_link(play, loop, time_slider, speedup)

        proxy.visible = None
        if self._system is None:
            resim.visible = None

        def on_resim_click(button):
            play.value = False
            self._resimulate_system()
            self._update_system_widgets()
        resim.on_click(on_resim_click)

    def _generate_pythreejs_meshes(self):
        self._meshes = []
        for obj in self.visualization_frames:
            geometry = _create_pythreejs_geometery(obj.shape)
            # TODO: set material
            material = p3js.LambertMaterial(color=obj.shape.color)
            mesh = p3js.Mesh(geometry=geometry, material=material)
            mesh.visualization_frame = obj
            self._meshes.append(mesh)

    def _generate_mesh_trajectories(self):
        for i, mesh in enumerate(self._meshes):
            position = []
            quaternion = []
            for transform in mesh.visualization_frame._visualization_matrix:
                position.append(transform[12:15])
                R = transform[0:3] + transform[4:7] + transform[8:11]
                quaternion.append(quaternion_from_rotation(R))
            self._trajectory_links[i].position = position
            self._trajectory_links[i].quaternion = quaternion

    def _update_system_widgets(self):
        # this should only be called if the system has been resimulated
        # which requires the system attribute to be set

        # setting the constants widget properties shouldn't be necessary
        # as we use the values in the widget when we call resimulate
        times = self._system.times
        t0 = times[0]
        tf = times[-1]
        dt = (tf - t0) / (len(times) - 1)

        # time step may not match as we use np.linspace to calculate the time
        # array
        self._time_widget.children[2].value = dt

        time_slider = self._play_widget.children[1]
        time_slider.min = t0
        time_slider.max = tf
        time_slider.step = dt
        time_slider.value = t0

    def _resimulate_system(self):
        # get parameter changes
        for w in self._constants_widget.children:
            self._system.constants[w._symbol] = w.value
        t = self._time_widget.children
        n = (t[1].value - t[0].value) / (t[2].value) + 1
        self._system.times = np.linspace(t[0].value, t[1].value, n)
        self._generate_simulation_dict()
        self._generate_mesh_trajectories()

    def display_pythreejs(self):
        if not hasattr(self, 'simulation_info'):
            self._generate_simulation_dict()

        self._generate_pythreejs_meshes()
        self._generate_time_widget()
        self._generate_constants_widget()
        self._generate_play_widget()

        self._trajectory_links = []
        for mesh in self._meshes:
            self._trajectory_links.append(
                    trajectory_link(self._play_widget.children[1], mesh))
        self._generate_mesh_trajectories()

        children = self._meshes + [p3js.AmbientLight(color=0x777777)]
        # TODO: lights
        # TODO: cameras
        c = p3js.PerspectiveCamera(
                position=[0, 15, 15], up=[0, 0, 1],
                children=[p3js.DirectionalLight(color='white',
                                                position=[3, 5, 1],
                                                intensity=0.5)])
        self._renderer = p3js.Renderer(
                camera=c,
                scene=p3js.Scene(children=children),
                controls=[p3js.OrbitControls(controlling=c)],
                width=700,
                height=400)
        self._widget = widgets.VBox()
        self._widget.children = (self._constants_widget,
                                 self._time_widget,
                                 self._play_widget,
                                 self._renderer)
        display(self._widget)


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
        The optional keyword arguments are the same as those in the
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

    def create_static_html(self, overwrite=False, silent=False, prefix=None):
        """Creates a directory named ``pydy-visualization`` in the current
        working directory which contains all of the HTML, CSS, Javascript,
        and json files required to run the vizualization application. To run
        the application, navigate into the ``pydy-visualization`` directory
        and start a webserver from that directory, e.g.::

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
            If True, the directory named ``pydy-visualization`` in the
            current working directory will be overwritten.
        silent : boolean, optional, default=False
            If True, no messages will be displayed to STDOUT.
        prefix : string, optional
            An optional prefix for the json data files.

        """

        pydy_dir = os.path.join(os.getcwd(), self.pydy_directory)

        if os.path.exists(pydy_dir) and overwrite is False:
            msg = ("The '{}' directory already exists. Would "
                   "you like to overwrite the contents? [y|n]\n")
            ans = raw_input(msg.format(self.pydy_directory))
            if ans == 'y':
                overwrite = True
            else:
                if not silent:
                    print("Aborted!")
                return

        # Copy all of the HTML/CSS/JS files from the source tree into the
        # local directory.
        if not silent:
            print("Copying static data.")
        src = os.path.join(os.path.dirname(__file__), 'static')
        distutils.dir_util.copy_tree(src, pydy_dir)

        # Add the two json files to the directory.
        if not silent:
            print("Copying Simulation data.")
        self._generate_json(directory=pydy_dir, prefix=prefix)

        if not silent:
            msg = ("To view the visualzation, run `python -m "
                   "SimpleHTTPServer` from the `{}` directory.")
            print(msg.format(self.pydy_directory))

    def remove_static_html(self, force=False):
        """Removes the ``static`` directory from the current working
        directory.

        Parameters
        ----------
        force : boolean, optional, default=False
            If true, no warning is issued before the removal of the
            directory.

        """

        pydy_dir = os.path.join(os.getcwd(), self.pydy_directory)

        if not os.path.exists(pydy_dir):
            print("All Done!")
            return

        if not force:
            msg = ("Are you sure you would like to delete the '{}' "
                   "directory? [y|n]\n")
            ans = raw_input(msg.format(self.pydy_directory))
            if ans == 'y':
                force = True

        if force:
            distutils.dir_util.remove_tree(pydy_dir)
            print("All Done!")
        else:
            print("Aborted!")

    def display(self):
        """Displays the scene in the default web browser."""
        self.create_static_html(overwrite=True, silent=True)
        resource_dir = os.path.join(os.getcwd(), self.pydy_directory)
        server = Server(scene_file=self._scene_json_file,
                        directory=resource_dir)
        server.run_server()

    def _rerun_button_callback(self, btn):
        """Callback for the "Rerun Simulation" button. When executed the
        parameter values are collected from the text input widgets and used
        in a new simulation of the model."""

        btn._dom_classes = ['btn-info', 'active', 'disabled']

        btn.description = 'Rerunning Simulation...'

        original_scene_file = self._scene_json_file
        original_sim_file = self._simulation_json_file
        original_constants = self.system.constants.copy()
        try:
            self.system.constants = {s: w.value for s, w in
                                     self._constants_text_widgets.items()}
            pydy_dir = os.path.join(os.getcwd(), self.pydy_directory)
            self._generate_json(directory=pydy_dir)
        except:
            print('Simulation rerun failed, using previous simulation data.')
            # If the simulation fails for any reason we revert everything
            # back to the previous state, including filling the text widgets
            # with the previous values of the constants. The _scene_info and
            # _simulation_data dicts may be in a bad state, but that should
            # be ok, because generate_visualiation_json will have to be run
            # again for anything new to happen.
            self._scene_json_file = original_scene_file
            self._simulation_json_file = original_sim_file
            self.system.constants = original_constants
            self._fill_constants_widgets()

        js_tmp = 'jQuery("#json-input").val("{}/{}");'
        js = js_tmp.format(self.pydy_directory, self._scene_json_file)
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

        self.create_static_html(overwrite=True, silent=True)

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

        ipython_static_url = os.path.relpath(self.pydy_directory, os.getcwd())
        ip_html_file = os.path.join(ipython_static_url, "index_ipython.html")
        with open(ip_html_file, 'r') as html_file:
            html = html_file.read()

        html = html.format(static_url=ipython_static_url,
                           load_url=os.path.join(ipython_static_url,
                                                 self._scene_json_file))

        self._html_widget = widgets.HTML(value=html)

        display(self._html_widget)
