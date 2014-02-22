====
PyDy
====

PyDy_, short for Python Dynamics, is a tool kit written in and accessed by the
Python programming language that utilizes an array of scientific tools to study
multibody dynamics. The goal is to have a modular framework which utilizes a
variety of tools that can provide the user with their desired workflow,
including:

- Model construction
- Equation of motion generation
- Simulation
- Visualization
- Publication

We started by building the SymPy_ `mechanics package`_ which provides an API
for building models and generating the symbolic equations of motion for complex
multibody systems and have more recently developed two packages, pydy-code-gen_
and pydy-viz_, for simulation and visualization of the models. The remaining
tools currently used in the PyDy workflow are popular scientific Python
packages such as NumPy_, SciPy_, IPython_, and matplotlib_ (i.e. the SciPy
stack) which provide additional code for numerical analyses, simulation, and
visualization.

Installation
============

The PyDy workflow generally depends on these Python packages:

SciPy Stack

- SymPy_ >= 0.7.2
- NumPy_ >= 1.6.1
- SciPy_ >= 0.9.0
- matplotlib_ >= 0.99.0
- IPython_ >= 0.13.0

PyDy Stack

- pydy-code-gen_ >= 0.1.0
- pydy-viz_ >= 0.1.0

It's best to install the dependencies from the SciPy Stack using the
instructions_ provided on the SciPy website.

Once you have all of the SciPy Stack dependencies you can simply install the
PyDy Stack with pip::

   $ pip install pydy

Or download the source and run::

   $ python setup.py install

For system wide installs you will need root permissions (perhaps prepend
commands with ``sudo``).

Note that the PyDy package is currently a simple wrapper to pydy-code-gen and
pydy-viz that provides a common namespace ``pydy``. These packages will likely
be merged into this package soon.

Usage
=====

Simply import the modules and functions when in a Python interpreter::

   >>> from sympy import symbols
   >>> from sympy.physics import mechanics
   >>> from pydy import codegen, visualization

Related Packages
================

- https://github.com/cdsousa/sympybotics
- https://pypi.python.org/pypi/Hamilton
- https://pypi.python.org/pypi/arboris
- https://pypi.python.org/pypi/PyODE
- https://pypi.python.org/pypi/odeViz
- https://pypi.python.org/pypi/ARS
- https://pypi.python.org/pypi/pymunk

.. _PyDy: http://pydy.org
.. _SymPy: http://sympy.org
.. _mechanics package: http://docs.sympy.org/latest/modules/physics/mechanics/index.html
.. _NumPy: http://numpy.scipy.org
.. _SciPy: http://www.scipy.org/scipylib/index.html
.. _matplotlib: http://matplotlib.org
.. _IPython: http://ipython.org
.. _pydy-code-gen: https://pypi.python.org/pypi/pydy-code-gen
.. _pydy-viz: https://pypi.python.org/pypi/pydy-viz
.. _instructions: http://www.scipy.org/install.html
