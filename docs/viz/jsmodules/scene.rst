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

addLights
^^^^^^^^^

*args*: None

Adds the Lights 
loaded from the JSON file
onto the scene. 

_addIndividualObject
^^^^^^^^^^^^^^^^^^^^

*args*: JS object, { object }

Adds a single geometry object
which is taken as an argument
to this function.

_addIndividualCamera
^^^^^^^^^^^^^^^^^^^^

*args*: JS object, { object }

Adds a single camera object
which is taken as an argument
to this function.

_addIndividualLight
^^^^^^^^^^^^^^^^^^^^

*args*: JS object, { object }

Adds a single light object
which is taken as an argument
to this function.

runAnimation
^^^^^^^^^^^^

*args*: None

This function iterates over the
the simulation data to render them
on the canvas.

setAnimationTime
^^^^^^^^^^^^^^^^

*args*: time, (float)

Takes a time value as the argument
and renders the simulation data 
corresponding to that time value.

stopAnimation
^^^^^^^^^^^^^

*args*: None

Stops the animation, and
sets the current time value to initial.

_removeAll
^^^^^^^^^^

*args*: None

Removes all the geometry elements
added to the scene from the loaded scene
JSON file. Keeps the default elements, i.e.
default axis, camera and light.

_blink
^^^^^^

*args*: id, (int)
Blinks the geometry element.
takes the element simulation_id as the 
argument and blinks it until some event is
triggered(UI button press)
