#!/home/tarun/anaconda/bin/python
from numpy import pi

from pydy.viz.shapes import Circle
from pydy.viz.scene import Scene
from pydy.viz.visualization_frame import VisualizationFrame

from simulate import *

# Create geometry
# ===============

# Each link in the pendulum is visualized with a cylinder, and a sphere at its
# far end.

disc = Circle(name='disc', radius=r, color="red")

disc_viz_frame = VisualizationFrame("disc_frame", R, Dmc, disc)

scene = Scene(N, C, disc_viz_frame)


# Create the visualization
# ========================

scene.generate_visualization_json(coordinates + speeds, constants.keys(), x,
        constants.values()) 
print len(x)
scene.create_static_html()
