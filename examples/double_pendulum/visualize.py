"""
This file will use pydy.viz to visualize the double pendulum.  Run this script
via a command line:

    $ python visualization.py
    
Note that you must add the pydy/viz directory to your PYTHONPATH for the
visualization to work:

    $ export PYTHONPATH=$PYTHONPATH:/path/to/pydy/viz

On my machine, I use:

    $ export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python2.7/dist-packages/pydy-0.1.0-py2.7.egg/pydy/viz

"""

from numpy import pi

from pydy.viz.shapes import Cylinder, Sphere
from pydy.viz.scene import Scene
from pydy.viz.visualization_frame import VisualizationFrame

from simulate import *


# Create geometry
# ===============

# Each link in the pendulum is visualized with a cylinder, and a sphere at its
# far end.
link = Cylinder(name='link', radius=0.5, length=constants[l], color='red')
sphere = Sphere(name='sphere', radius=1.0)

# TODO notes.
linkP_viz_frame = VisualizationFrame('linkP',
        A.orientnew('frameP', 'Axis', [0.5 * pi, N.z]),
        O.locatenew('originP', 0.5 * l * A.x),
        link)
# TODO notes.
linkR_viz_frame = VisualizationFrame('linkR',
        B.orientnew('frameR', 'Axis', [0.5 * pi, N.z]),
        P.locatenew('originP', 0.5 * l * B.x),
        link)

sphereP_viz_frame = VisualizationFrame('sphereP', N, P, sphere)
sphereR_viz_frame = VisualizationFrame('sphereR', N, R, sphere)


# Construct the scene
# ===================

# TODO world frame...
scene = Scene(N.orientnew('world', 'Axis', [0.5 * pi, N.z]), O,
        linkP_viz_frame, linkR_viz_frame, sphereP_viz_frame, sphereR_viz_frame)


# Create the visualization
# ========================

scene.generate_visualization_json(coordinates + speeds, constants.keys(), x,
        constants.values())

scene.display()
