These scripts generate functions to compute the lateral tire-ground constraint
forces in the Carvallo-Whipple model. This assumes no-slip point contact. The
coordinate and parameter definitions follow those in [Moore2012]_.

.. [Moore2012] Moore, Jason K. "Human Control of a Bicycle." Doctor of
   Philosophy, University of California, 2012.
   http://moorepants.github.io/dissertation.

First run the Python file::

   $ python whipple.py

This generates 4 Octave files:

- ``eval_holonomic.m``
- ``eval_dep_speeds.m``
- ``eval_dep_speeds_derivs.m``
- ``eval_lat_forces.m``

Once those files are available, the ``lateral_tire_forces.m`` provides a
function to compute the forces given the kinematics of the system. Run::

   $ octave try_lateral_forces.m

To try the function with some sample values.
