JavaScript functions Reference
-------------------------------

Canvas
======
Canvas is the base class for handling all the animation and 
visualization generation.

canvas/initialize.js
====================
Constructor:
^^^^^^^^^^^^

This function acts as a class constructor for Canvas class
It takes the JSON Object variable as the argument, which contains
all the data in the JSON format.
It binds onClick methods of certain Divs on the frontend
with some Canvas.prototype functions.

Canvas.prototype.initialize
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function initializes the starting canvas, on which
all the visualizations are drawn. 

It adds following to the canvas:
  - A Primary Camera
  - Primary Trackball Controls
  - A Primary Light
  - Axes
  - Grid
  - A Div for displaying total number of frames.
  - A Div for displaying the current frame animation is
    running on.

canvas/addObjects.js
====================
Canvas.prototype.addControls
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function initializes the Primary Controls,
which were defined in Canvas.prototype.initialize function.

It generates a controlsID, which contains the return value 
of requestAnimationFrame, and can be used to call 
cancelAnimationFrame, for stopping mouse controlled animation.

Canvas.prototype.resetControls
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function simply calls the controls reset method
for the Primary Controls(canvas.prototype.primaryControls).

Canvas.prototype.addCameras
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function parses the JSON Object for cameras
and adds them to the scene.
All the cameras are stored in a Canvas.cameras object,
which is an instance of THREE.Object3D();

Canvas.prototype.addLights
^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function parses the JSON Object for lights
and adds them to the scene.
All the lights are stored in a Canvas.lights object,
which is an instance of THREE.Object3D();

Canvas.prototype.addFrames
^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function parses the JSON Object for frames
and adds them to the scene.
All the frames are stored in a Canvas.frames object,
which is an instance of THREE.Object3D();

canvas/animate.js
=================
Canvas.prototype.startAnimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function kick starts the animation.
It iterates over the frames and apply transformation matrices
from Simulation Matrix of that frame, iteratively.
by default animation is done for a single loop,
which can be changed to looped by the check button from the UI.

Canvas.prototype.pauseAnimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function pauses the animation, but retains the
current animation frame.

Canvas.prototype.stopAnimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This prototype function stops the animation, and resets 
current animation frame to 0.
