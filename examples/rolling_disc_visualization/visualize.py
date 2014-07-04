from numpy import pi

from pydy.viz.shapes import Cylinder
from pydy.viz.scene import Scene
from pydy.viz.visualization_frame import VisualizationFrame

from simulate import *

# Create geometry
# ===============

disc = Cylinder(name='disc', length=0.01, radius=r, color="red")

disc_viz_frame = VisualizationFrame("disc_frame", R, Dmc, disc)

# In the derivation, the "ground" is the xy plane. However, the ground plane in
# the visualizer is the xz plane. We rotate everything 90 degrees around the
# Newtonian x axis. so that the disc starts standing on the "ground" in the
# visualizer.
world_frame = N.orientnew('world', 'Axis', [0.5 * pi, N.x])
scene = Scene(world_frame, No, disc_viz_frame)


# Create the visualization
# ========================

scene.generate_visualization_json(coordinates + speeds, constants.keys(), x,
        constants.values()) 
print len(x)
scene.create_static_html()
