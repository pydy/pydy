====
PyDy
====

|pypi| |binstar| |rtd-docs| |travis-build| |Gitter Badge|

.. |pypi| image:: https://img.shields.io/pypi/v/pydy.svg
   :target: https://pypi.python.org/pypi/pydy
   :alt: Latest Released Version

.. |binstar| image:: https://binstar.org/pydy/pydy/badges/version.svg
   :target: https://binstar.org/pydy/pydy

.. |rtd-docs| image:: https://readthedocs.org/projects/pydy/badge/?version=latest
   :target: http://pydy.readthedocs.org/en/latest
   :alt: Latest documentation

.. |travis-build| image:: https://travis-ci.org/pydy/pydy.png?branch=master
   :target: https://travis-ci.org/pydy/pydy

.. |Gitter Badge| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/pydy/pydy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

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
multibody systems. More recently we developed two packages, `pydy.codegen`
and `pydy.viz`, for simulation and visualization of the models, respectively.
This Python package contains these two packages and other tools for working
with SymPy mechanics. The remaining tools currently used in the PyDy workflow
are popular scientific Python packages such as NumPy_, SciPy_, IPython_, and
matplotlib_ (i.e. the SciPy stack) which provide additional code for numerical
analyses, simulation, and visualization.

.. _SymPy: http://sympy.org
.. _mechanics package: http://docs.sympy.org/latest/modules/physics/mechanics/index.html
.. _NumPy: http://numpy.scipy.org
.. _SciPy: http://www.scipy.org/scipylib/index.html
.. _IPython: http://ipython.org
.. _matplotlib: http://matplotlib.org

Installation
============

PyDy has hard dependencies on the following software\ [#]_:

- Python >= 2.7, >= 3.3
- setuptools\ [#]_
- NumPy_ >= 1.7 (>= 1.9 with Python >= 3.3)
- SciPy_ >= 0.11.0 (>= 0.14.0 with Python >= 3.3)
- SymPy_ >= 0.7.4.1 (>= 0.7.5 with Python >= 3.3)

.. [#] We only test PyDy with these minimum dependencies. Previous versions may
   work.
.. [#] setuptools >= 8.0 is required if development versions of SymPy are used.

PyDy has optional dependencies on these packages:

- IPython_ >= 0.3.0
- Theano_ >= 0.6.0 (>= 0.7.0 with Python >= 3.3)
- Cython_ >= 0.17 (>= 0.20.1 with Python >= 3.3)

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

Once the dependencies are installed, the latest stable version of the package
can be downloaded from PyPi::

   $ wget https://pypi.python.org/packages/source/p/pydy/pydy-0.2.1.tar.gz

and extracted and installed\ [#]_::

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
                specifieds={force: lambda x, t: sin(t)},
                initial_conditions={position: 0.1, speed: -1.0},
                times=linspace(0.0, 10.0, 1000))

Integrate the equations of motion to get the state trajectories::

   y = sys.integrate()

Plot the results::

   import matplotlib.pyplot as plt

   plt.plot(sys.times, y)
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

For the JavaScript tests the Jasmine and blanket.js libraries are used. Both
of these libraries are included in pydy.viz with the source. To run the
JavaScript tests::

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

Questions, Bugs, Feature Requests
=================================

If you have any question about installation, usage, etc, feel free send a
message to our public `mailing list`_ or visit our `Gitter chatroom`_.

.. _mailing list: http://groups.google.com/group/pydy
.. _Gitter chatroom: https://gitter.im/pydy/pydy

If you think there’s a bug or you would like to request a feature, please open
an `issue`_ on Github.

.. _issue: https://github.com/pydy/pydy/issues

Release Notes
=============

0.3.0
-----

User Facing
~~~~~~~~~~~

- Introduced conda builds and binstar support. [PR `#219`_]
- Dropped support for IPython < 3.0. [PR `#237`_]
- Added support Python 3.3 and 3.4. [PR `#229`_]
- Bumped up the minimum dependencies for NumPy, SciPy, and Cython [PR `#233`_].
- Removed the partial implementation of the Mesh shape. [PR `#172`_]
- Overhauled the code generation package to make the generators more easily
  extensible and to improve simulation speed. [PR `#113`_]
- The visualizer has been overhauled as part of Tarun Gaba's 2014 GSoC
  internship [PR `#82`_]. Here are some of the changes:

  - The JavaScript is now handled by AJAX and requires a simple server.
  - The JavaScript has been overhauled and now uses prototype.js for object
    oriented design.
  - The visualizer can now be loaded in an IPython notebook via IPython's
    widgets using ``Scene.display_ipython()``.
  - A slider was added to manually control the frame playback.
  - The visualization shapes' attributes can be manipulated via the GUI.
  - The scene json file can be edited and downloaded from the GUI.
  - pydy.viz generates two JSONs now (instead of one in earlier versions). The
    JSON generated from earlier versions will **not** work in the new version.
  - Shapes can now have a material attribute.
  - Model constants can be modified and the simulations can be rerun all via
    the GUI.
  - Switched from socket based server to python's core SimpleHTTPServer.
  - The server has a proper shutdown response [PR `#241`_]

- Added a new experimental System class and module to more seamlessly manage
  integrating the equations of motion. [PR `#81`_]

.. _#241: https://github.com/pydy/pydy/pull/241
.. _#237: https://github.com/pydy/pydy/pull/237
.. _#229: https://github.com/pydy/pydy/pull/229
.. _#233: https://github.com/pydy/pydy/pull/233
.. _#219: https://github.com/pydy/pydy/pull/219
.. _#172: https://github.com/pydy/pydy/pull/172
.. _#113: https://github.com/pydy/pydy/pull/113
.. _#82: https://github.com/pydy/pydy/pull/82
.. _#81: https://github.com/pydy/pydy/pull/81

Development
~~~~~~~~~~~

- Switched to a conda based Travis testing setup. [PR `#231`_]
- When using older SymPy development versions with non-PEP440 compliant version
  identifiers, setuptools < 8 is required. [PR `#166`_]
- Development version numbers are now PEP 440 compliant. [PR `#141`_]
- Introduced pull request checklists and CONTRIBUTING file. [PR `#146`_]
- Introduced light code linting into Travis. [PR `#148`_]

.. _#231: https://github.com/pydy/pydy/pull/231
.. _#166: https://github.com/pydy/pydy/pull/166
.. _#141: https://github.com/pydy/pydy/pull/141
.. _#146: https://github.com/pydy/pydy/pull/146
.. _#148: https://github.com/pydy/pydy/pull/148

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
