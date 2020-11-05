These scripts generate functions to compute the lateral tire-ground constraint
forces in the Carvallo-Whipple bicycle model. This assumes no-slip point
contact. The coordinate and parameter definitions follow those in [Moore2012]_.

.. [Moore2012] Moore, Jason K. "Human Control of a Bicycle." Doctor of
   Philosophy, University of California, 2012.
   http://moorepants.github.io/dissertation.

Dependencies
============

The Python file needs SymPy and PyDy. If you install Anaconda or Miniconda you
can install the latest versions of these with::

   $ conda install -c conda-forge sympy pydy

If you use pip for installing packages you can use::

   $ pip install sympy pydy

Octave or Matlab is needed for the ``m`` files.

Usage
=====

First run the Python file::

   $ python whipple.py

This generates 4 Octave files:

- ``eval_holonomic.m``
- ``eval_dep_speeds.m``
- ``eval_dep_speeds_derivs.m``
- ``eval_lat_forces.m``

These 4 files are also checked into the Git repository and you can skip
regenerating them with Python.

Once those files are available, the ``lateral_tire_forces.m`` provides a
function to compute the forces given the kinematics of the system. To execute
that function with some sample numerical values run::

   $ octave try_lateral_forces.m
