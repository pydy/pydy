**Warning: This package is in development and not ready for use by users. We
are migrating to a simpy one package installer for the PyDy software suite.***

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
an API for building models and generating the equations of motion for complex
multibody systems. The remaining tools used in the PyDy workflow are popular
scientific Python packages such as NumPy, SciPy, IPython, and matplotlib which
provide code for numerical analyses, simulation, and visualization.

Installation
============

Get the requirements:

- SymPy >= 0.7.2 (best to use the dev version as the development is active)
- NumPy
- SciPy
