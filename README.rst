====
PyDy
====

.. image:: https://img.shields.io/pypi/dm/pydy.svg
    :target: https://pypi.python.org/pypi/pydy/
    :alt: Latest Version

.. image:: https://travis-ci.org/pydy/pydy.png?branch=master
   :target: https://travis-ci.org/pydy/pydy

PyDy_, short for Python Dynamics, is a tool kit written in the Python
programming language that utilizes an array of scientific programs to enable
the study of multibody dynamics. The goal is to have a modular framework and
eventually a physics abstraction layer which utilizes a variety of backends
that can provide the user with their desired workflow, including:

.. _PyDy: http://pydy.org

- Model specification
- Equation of motion generation
- Simulation
- Visualization
- Publication

We started by building the SymPy_ `mechanics package`_ which provides an API
for building models and generating the symbolic equations of motion for complex
multibody systems and have more recently developed two packages, `pydy.codegen`
and `pydy.viz`, for simulation and visualization of the models, respectively.
The remaining tools currently used in the PyDy workflow are popular scientific
Python packages such as NumPy_, SciPy_, IPython_, and matplotlib_ (i.e. the
SciPy stack) which provide additional code for numerical analyses, simulation,
and visualization.

.. _SymPy: http://sympy.org
.. _mechanics package: http://docs.sympy.org/latest/modules/physics/mechanics/index.html
.. _NumPy: http://numpy.scipy.org
.. _SciPy: http://www.scipy.org/scipylib/index.html
.. _IPython: http://ipython.org
.. _matplotlib: http://matplotlib.org

Installation
============

PyDy has hard dependencies on the following software:

- Python >= 2.7, < 3
- setuptools
- SymPy_ >= 0.7.4.1
- NumPy_ >= 1.6.1
- SciPy_ >= 0.9.0
- IPython_ >= 0.13.0

PyDy has optional dependencies on these packages:

- Cython_ >= 0.15.1
- Theano_ >= 0.6.0

.. _Theano: http://deeplearning.net/software/theano/
.. _Cython: http://cython.org/

The examples may require these dependencies:

- matplotlib_

It's best to install the SciPy Stack dependencies using the instructions_
provided on the SciPy website. We recommend the conda_ package manager and the
Anaconda_ distribution for easy cross platform installation.

.. _instructions: http://www.scipy.org/install.html
.. _conda: http://conda.pydata.org/
.. _Anaconda: http://docs.continuum.io/anaconda/

Once the dependencies are installed, the package can be downloaed from PyPi::

   $ wget https://pypi.python.org/packages/source/p/pydy/pydy-0.2.1.tar.gz

and extracted and installed\ [#]::

   $ tar -zxvf pydy-0.2.1.tar.gz
   $ cd pydy-0.2.1
   $ python setup.py install

Or if you have the pip package manager installed you can simply type::

   $ pip install pydy

.. [#] For system wide installs you may need root permissions (perhaps prepend
   commands with ``sudo``).

Usage
=====

This is an example of a simple one degree of freedom system: a mass under the
influence of a spring, damper, gravity and an external force::


   / / / / / / / / /
   -----------------
     |    |     |   | g
     \   | |    |   V
   k /   --- c  |
     |    |     | x, v
    --------    V
    |  m   | -----
    --------
       | F
       V

Derive the system::

   from sympy import symbols
   import sympy.physics.mechanics as me

   mass, stiffness, damping, gravity = symbols('m, k, c, g')

   position, speed = me.dynamicsymbols('x v')
   positiond = me.dynamicsymbols('x', 1)
   force = me.dynamicsymbols('F')

   ceiling = me.ReferenceFrame('N')

   origin = me.Point('origin')
   origin.set_vel(ceiling, 0)

   center = origin.locatenew('center', position * ceiling.x)
   center.set_vel(ceiling, speed * ceiling.x)

   block = me.Particle('block', center, mass)

   kinematic_equations = [speed - positiond]

   force_magnitude = mass * gravity - stiffness * position - damping * speed + force
   forces = [(center, force_magnitude * ceiling.x)]

   particles = [block]

   kane = me.KanesMethod(ceiling, q_ind=[position], u_ind=[speed],
                         kd_eqs=kinematic_equations)
   kane.kanes_equations(forces, particles)

Create a system to manage integration and specify numerical values for the
constants and specified quantities. Here, we specify sinusoidal forcing::

   from numpy import array, linspace, sin
   from pydy.system import System

   sys = System(kane,
                constants={mass: 1.0, stiffness: 1.0,
                           damping: 0.2, gravity: 9.8},
                specified=lambda x, t: sin(t),
                initial_conditions=array([0.1, -1.0]))

Now generate the function needed for numerical evaluation of the ODEs::

   sys.generate_ode_function()

Integrate the equations of motion under the influence of a specified sinusoidal
force::

   t = linspace(0.0, 10.0, 1000)
   y = sys.integrate(t)

Plot the results::

   import matplotlib.pyplot as plt

   plt.plot(t, y)
   plt.legend((str(position), str(speed)))
   plt.show()

Documentation
=============

The documentation is hosted at http://pydy.readthedocs.org but you can also
build them from source using the following instructions.

To build the documentation you must install the dependencies:

- Sphinx_
- numpydoc_

.. _Sphinx: http://sphinx-doc.org/
.. _numpydoc: https://pypi.python.org/pypi/numpydoc

To build the HTML docs, run Make from within the ``docs`` directory::

   $ cd docs
   $ make html

You can then view the documentation from your preferred web browser, for
example::

   $ firefox _build/html/index.html

Modules and Packages
====================

Code Generation (codegen)
-------------------------

This package provides code generation facilities. It generates functions that
can numerically evaluate the right hand side of the ordinary differential
equations generated with sympy.physics.mechanics_ with three different
backends: SymPy's lambdify_, Theano, and Cython.

.. _sympy.physics.mechanics: http://docs.sympy.org/latest/modules/physics/mechanics
.. _lambdify: http://docs.sympy.org/latest/modules/utilities/lambdify.html#sympy.utilities.lambdify.lambdify

Models (models.py)
------------------

The models module provides some canned models of classic systems.

Systems (system.py)
-------------------

The System module provides a ``System`` class to manage simulation of a single
system.

Visualization (viz)
-------------------

This package provides tools to create 3D animated visualizations of the
systems. The visualizations utilize WebGL and run in a web browser. They can
also be embedded into an IPython notebook for added interactivity.

Development Environment
=======================

The source code is managed with the Git version control system. To get the
latest development version and access to the full repository, clone the
repository from Github with::

   $ git clone https://github.com/pydy/pydy.git

You should then install the dependencies for running the tests:

- nose_: 1.3.0
- phantomjs_: 1.9.0

.. _nose: https://nose.readthedocs.org
.. _phantomjs: http://phantomjs.org

Isolated Environments
---------------------

It is typically advantageous to setup a virtual environment to isolate the
development code from other versions on your system. There are two popular
environment managers that work well with Python packages: virtualenv and
conda_.

The following installation assumes you have virtualenvwrapper_ in addition to
virtualenv and all the dependencies needed to build the various packages::

   $ mkvirtualenv pydy-dev
   (pydy-dev)$ pip install numpy scipy cython nose theano sympy ipython[all]
   (pydy-dev)$ pip install matplotlib # make sure to do this after numpy
   (pydy-dev)$ git clone git@github.com:pydy/pydy.git
   (pydy-dev)$ cd pydy
   (pydy-dev)$ python setup.py develop

.. _virtualenvwrapper: https://pypi.python.org/pypi/virtualenvwrappe://pypi.python.org/pypi/virtualenvwrapper

Or with conda_::

   $ conda create -n pydy-dev setuptools numpy scipy ipython ipython-notebook cython nose theano sympy matplotlib
   $ source activate pydy-dev
   (pydy-dev)$ git clone git@github.com:pydy/pydy.git
   (pydy-dev)$ cd pydy
   (pydy-dev)$ python setup.py develop

The full Python test suite can be run with::

   (pydy-dev)$ nosetests

For the Javascript tests the Jasmine and blanket.js libraries are used. Both
of these libraries are included in pydy.viz with the source. To run the
Javascript tests::

   cd pydy/viz/static/js/tests && phantomjs run-jasmine.js SpecRunner.html && cd ../../../../../

Benchmark
=========

Run the benchmark to test the n-link pendulum problem with the various backends::

   $ python bin/benchmark_pydy_code_gen.py <max # of links> <# of time steps>

Related Packages
================

These are various related and similar Python packages:

- https://github.com/cdsousa/sympybotics
- https://pypi.python.org/pypi/Hamilton
- https://pypi.python.org/pypi/arboris
- https://pypi.python.org/pypi/PyODE
- https://pypi.python.org/pypi/odeViz
- https://pypi.python.org/pypi/ARS
- https://pypi.python.org/pypi/pymunk

Citation
========

If you make use of PyDy in your work or research, please cite us in your
publications or on the web. This citation can be used:

   Gilbert Gede, Dale L Peterson, Angadh S Nanjangud, Jason K Moore, and Mont
   Hubbard, "Constrained Multibody Dynamics With Python: From Symbolic Equation
   Generation to Publication", ASME 2013 International Design Engineering
   Technical Conferences and Computers and Information in Engineering
   Conference, 2013, `10.1115/DETC2013-13470
   <http://dx.doi.org/10.1115/DETC2013-13470>`_.

Release Notes
=============

0.3.0
-----

- Overhauled the code generation package to make the generators more easily
  extensible and to improve simluation speed.
- Added a new System class and module to more seamlessly manage integrating the
  equations of motion.

0.2.1
-----

- Unbundled unnecessary files from tar ball.

0.2.0
-----

- Merged pydy_viz, pydy_code_gen, and pydy_examples into the source tree.
- Added a method to output "static" visualizations from a Scene object.
- Dropped the matplotlib dependency and now only three.js colors are valid.
- Added joint torques to the n_pendulum model.
- Added basic examples for codegen and viz.
- Graceful fail if theano or cython are not present.
- Shapes can now use sympy symbols for geometric dimensions.
