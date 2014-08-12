DynamicsVisualizer
==================

DynamicsVisualizer is the main class for 
Dynamics Visualizer. It contains methods to 
set up a default UI, and maps buttons' 
`onClick` to functions.

_initialize
^^^^^^^^^^^
*args*: None

Checks whether the browser supports webGLs, and
initializes the DynamicVisualizer object.

isWebGLCompatible
^^^^^^^^^^^^^^^^^

*args*: None

Checks whether the browser used is
compatible for handling webGL based animations.
Requires external script: Modernizr.js

activateUIControls 
^^^^^^^^^^^^^^^^^^

*args*: None

This method adds functions to the UI buttons
It should be **strictly** called after the 
other DynamicsVisualizer sub-modules are loaded
in the browser, else certain functionality will 
be(not might be!) hindered.

loadUIElements
^^^^^^^^^^^^^^

*args*: None

This method loads UI elements 
which can be loaded only **after**
scene JSON is loaded onto canvas.

getBasePath
^^^^^^^^^^^

*args*: None

Returns the base path of the loaded Scene file.

getFileExtenstion
^^^^^^^^^^^^^^^^^

*args*: None

Returns the extension of
the uploaded Scene file.

getQueryString
^^^^^^^^^^^^^^

*args*: key

Returns the GET Parameter from url corresponding
to `key`