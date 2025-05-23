from numpy import pi
from pydy.viz.shapes import Cylinder, Plane
from pydy.viz.scene import Scene
from pydy.viz.visualization_frame import VisualizationFrame

from simulate import *

# Create geometry
# ===============

# The ground plane is centered at the origin and normal to the Z axis by
# default. It would be nice if the center of the ground plane could be at
# the center of the disc's path but that is trickier. You could create a new
# point in the ground plane with constants for the measure numbers then sub
# those in at this step.
ground_plane = Plane(color='blue',
                     length=2 * (max(x[:, 3]) - min(x[:, 3])),
                     width=2 * (max(x[:, 4]) - min(x[:, 4])))
ground_viz_frame = VisualizationFrame("ground", N, C, ground_plane)

disc = Cylinder(name='disc', radius=r, length=0.0001, color="red")
disc_viz_frame = VisualizationFrame("disc_frame", R, translating_com, disc)

scene = Scene(N, C, disc_viz_frame, ground_viz_frame)

# Create the visualization
# ========================

scene.generate_visualization_json(coordinates + speeds,
                                  constants.keys(),
                                  x,
                                  constants.values())
scene.display()
