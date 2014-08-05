DynamicsVisualizer.Scene
==============================

create
^^^^^^^^^^
*args*: None

This method creates the scene 
from the self.model
and renders it onto the canvas.

_createRenderer
^^^^^^^^^^^^^^^

*args*: None

Creates a webGL Renderer
 with a default background color.

_addDefaultLightsandCameras
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*args*: None

This method adds a default light
and a Perspective camera to the
initial visualization

_addAxes
^^^^^^^^

*args*: None

Adds a default system of axes
to the initial visualization.

_addTrackBallControls
^^^^^^^^^^^^^^^^^^^^^

*args*: None

Adds Mouse controls 
to the initial visualization
using Three's TrackballControls library.

_resetControls
^^^^^^^^^^^^^^

*args*: None

Resets the scene camera to 
the initial values(zoom, displacement etc.)

addObjects
^^^^^^^^^^

*args*: None

Adds the geometries 
loaded from the JSON file
onto the scene. The file is 
saved as an object in self.model
and then rendered to canvas with this
function.

addCameras
^^^^^^^^^^

*args*: None

Adds the cameras 
loaded from the JSON file
onto the scene. The cameras
can be switched during animation
from the `switch cameras` UI button.





