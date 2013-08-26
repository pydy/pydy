**Warning: This package is in development and not ready for use by end users.
We are migrating to a simple one package installer for the PyDy software
suite.**

====
PyDy
====

PyDy, short for Python Dynamics, is a tool kit written in and accessed by the
Python programming language that utlizes an array of scientific tools to study
multibody dynamics. The goal is to have a modular framework which utlizes a
variety of tools that can provide the user with their desired workflow,
including:

- Model construction
- Equation of motion generation
- Simulation
- Visualization
- Publication

Currently we've started by building the SymPy mechanics package which provides
an API for building models and generating the symbolic equations of motion for
complex multibody systems. The remaining tools currently used in the PyDy
workflow are popular scientific Python packages such as NumPy, SciPy, IPython,
and matplotlib which provide code for numerical analyses, simulation, and
visualization.

Installation
============

PyDy depends on:

- SymPy >= 0.7.2 (best to use the dev version as the development is active)
- pydy-viz >= 0.1.0
- NumPy
- SciPy
- matplotlib
- IPython

If you have all of the dependencies or all the dependencies needed to build the
packages above you can simply install with easy_install or pip::

   $ easy_install PyDy

or::

   $ pip install PyDy

Or download the source and run::

   $ python setup.py install

For system wide installs, prepend with `sudo`.

Usage
=====

Simply import the modules and functions when in a Python interpreter::

   >>> from sympy import symbols
   >>> from pydy import mechanics, visualization
