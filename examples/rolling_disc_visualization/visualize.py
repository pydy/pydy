#!/home/tarun/anaconda/bin/python
from numpy import pi

from pydy.viz.shapes import Cylinder, Plane
from pydy.viz.scene import Scene
from pydy.viz.visualization_frame import VisualizationFrame

from simulate import *

# Create geometry
# ===============

# Each link in the pendulum is visualized with a cylinder, and a sphere at its
# far end.

disc = Cylinder(name='disc', radius=r, length=0.0001, color="red")

disc_viz_frame = VisualizationFrame("disc_frame", R, translating_com, disc)

ground_plane = Plane(color='blue', length=max(x[:, 3]), width=max(x[:, 4]))
ground_viz_frame = VisualizationFrame("ground", N, C, ground_plane)

scene = Scene(N, C, disc_viz_frame, ground_viz_frame)

# Create the visualization
# ========================

scene.generate_visualization_json(coordinates + speeds, constants.keys(), x,
        constants.values())
print len(x)
scene.create_static_html()
