DynamicsVisualizer.ParamEditor
==============================

openDialog
^^^^^^^^^^
*args*: id, (str)

This function takes object's id
as the argument, and populates the
edit objects dialog box.

applySceneInfo
^^^^^^^^^^^^^^

*args*: id, (str)

This object applies the changes made in
the edit objects dialog box to self.model
and then renders the model onto canvas.
It takes the id of the object as its argument.

_addGeometryFor
^^^^^^^^^^^^^^^

*args*: JS object,{ object }

Adds geometry info for a particular
object onto the edit objects dialog
box. Takes the object as the argument.

showModel
^^^^^^^^^

*args*: None

Updates the codemirror instance with
the updated model, and shows it in the UI.
