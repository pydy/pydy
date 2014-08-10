DynamicVisualizer.Parser
========================

loadScene
^^^^^^^^^

*args*: None

This method calls an ajax request on the 
JSON file and reads the scene info from 
the JSON file, and saves it as an object
at self.model.

loadSimulation
^^^^^^^^^^^^^^

*args*: None

This method loads the simulation data 
from the simulation JSON file. The data is
saved in the form of 4x4 matrices mapped to 
the simulation object id, at a particular time.

createTimeArray
^^^^^^^^^^^^^^^

*args*: None

Creates a time array from 
the information inferred from
simulation data.
