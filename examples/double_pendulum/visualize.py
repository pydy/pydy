"""
This file will use pydy.viz to visualize the double pendulum.  Run this script
via a command line:

    $ python visualization.py

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
link = Cylinder(name='link', radius=0.5, length=l, color='red')
sphere = Sphere(name='sphere', radius=1.0)

# By default, Cylinders are drawn so that their center is at the origin of the
# VisualizationFrame, and their axis is the y axis of the VisualizationFrame.
# We want the end of the Cylinder to be at the origin of the
# VisualizationFrame, and we want the Cylinder's axis to be aligned with the x
# axis of the VisualizationFrame. For these reasons, we must use the
# 'orientnew' and 'locatenew' methods to create new frames/points.
linkP_frame = A.orientnew('frameP', 'Axis', [0.5 * pi, N.z])
linkP_origin = O.locatenew('originP', 0.5 * l * A.x)
linkP_viz_frame = VisualizationFrame('linkP', linkP_frame, linkP_origin, link)

linkR_frame = B.orientnew('frameR', 'Axis', [0.5 * pi, N.z])
linkR_origin = P.locatenew('originP', 0.5 * l * B.x)
linkR_viz_frame = VisualizationFrame('linkR', linkR_frame, linkR_origin, link)

sphereP_viz_frame = VisualizationFrame('sphereP', N, P, sphere)
sphereR_viz_frame = VisualizationFrame('sphereR', N, R, sphere)


# Construct the scene
# ===================

# We want gravity to be directed downwards in the visualization. Gravity is in
# the -x direction. By default, the visualization uses the xz plane as the
# ground plane. Thus, gravity is contained in the ground plane. However, we
# want gravity to point in the -y direction in the visualization. To achieve
# this, we create a world frame that is rotated +90 degrees about the N frame's
# z direction.
world_frame = N.orientnew('world', 'Axis', [0.5 * pi, N.z])
scene = Scene(world_frame, O,
        linkP_viz_frame, linkR_viz_frame, sphereP_viz_frame, sphereR_viz_frame)


# Create the visualization
# ========================

scene.generate_visualization_json_system(sys)

if __name__ == "__main__":
    try: #If called from inside notebook,
        scene.display_ipython()
    except:#If called from interpreter
        scene.display()
