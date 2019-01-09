====
PyDy
====

|pypi| |anaconda| |rtd-docs| |travis-build| |appveyor| |gitter|

.. |pypi| image:: https://img.shields.io/pypi/v/pydy.svg
   :target: https://pypi.python.org/pypi/pydy
   :alt: Latest Released Version

.. |anaconda| image:: https://anaconda.org/pydy/pydy/badges/version.svg
   :target: https://anaconda.org/pydy/pydy

.. |rtd-docs| image:: https://readthedocs.org/projects/pydy/badge/?version=stable
   :target: https://pydy.readthedocs.io/en/stable/?badge=stable
   :alt: Documentation Status

.. |travis-build| image:: https://travis-ci.org/pydy/pydy.png?branch=master
   :target: https://travis-ci.org/pydy/pydy

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/orj87gyb0c1wqc6j/branch/master?svg=true
   :target: https://ci.appveyor.com/project/moorepants/pydy/branch/master

.. |gitter| image:: https://badges.gitter.im/Join%20Chat.svg
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
multibody systems. More recently we developed two packages, `pydy.codegen` and
`pydy.viz`, for simulation and visualization of the models, respectively.  This
Python package contains these two packages and other tools for working with
mathematical models generated from SymPy mechanics. The remaining tools
currently used in the PyDy workflow are popular scientific Python packages such
as NumPy_, SciPy_, IPython_, Jupyter_, ipywidgets_, and matplotlib_ (i.e. the
SciPy stack) which provide additional code for numerical analyses, simulation,
and visualization.

.. _SymPy: http://sympy.org
.. _mechanics package: http://docs.sympy.org/latest/modules/physics/mechanics/index.html
.. _NumPy: http://numpy.scipy.org
.. _SciPy: http://www.scipy.org/scipylib/index.html
.. _IPython: http://ipython.org
.. _Jupyter: http://jupyter.org
.. _ipywidgets: https://pypi.python.org/pypi/ipywidgets
.. _matplotlib: http://matplotlib.org

Installation
============

PyDy has hard dependencies on the following software\ [#]_:

.. [#] We only test PyDy with these minimum dependencies; these module versions
       are provided in the Ubuntu 16.04 packages. Previous versions may work.

- 2.7 <= Python < 3.0 or Python >= 3.5
- setuptools >= 20.7.0
- NumPy_ >= 1.11.0
- SciPy_ >= 0.17.1
- SymPy_ >= 0.7.6.1
- PyWin32 >= 219 (Windows Only)

PyDy has optional dependencies on these packages:

- 4.0.0 <= `Jupyter Notebook`_ < 5.0.0
- 4.0.0 <= ipywidgets_ < 5.0.0
- Theano_ >= 0.8.0
- Cython_ >= 0.23.4

.. _Theano: http://deeplearning.net/software/theano/
.. _Cython: http://cython.org/
.. _Jupyter Notebook: https://pypi.python.org/pypi/notebook

The examples may require these dependencies:

- matplotlib_ >= 1.5.1
- version_information_

.. _version_information: https://pypi.python.org/pypi/version_information

It's best to install the SciPy Stack dependencies using the instructions_
provided on the SciPy website first. We recommend the conda_ package manager
and the Anaconda_ distribution for easy cross platform installation.

.. _instructions: http://www.scipy.org/install.html
.. _conda: http://conda.pydata.org/
.. _Anaconda: http://docs.continuum.io/anaconda/

Once the dependencies are installed, the latest stable version of the package
can be downloaded from PyPi\ [#]_::

   $ wget https://pypi.python.org/packages/source/p/pydy/pydy-X.X.X.tar.gz

.. [#] Change ``X.X.X`` to the latest version number.

and extracted and installed\ [#]_::

   $ tar -zxvf pydy-X.X.X.tar.gz
   $ cd pydy-X.X.X
   $ python setup.py install

.. [#] For system wide installs you may need root permissions (perhaps prepend
   commands with ``sudo``).

Or if you have the pip package manager installed you can simply type::

   $ pip install pydy

Or if you have conda you can type::

   $ conda install -c conda-forge pydy

Also, a simple way to install all of the optional dependencies is to install
the ``pydy-examples`` metapackage using conda::

   $ conda install -c pydy pydy-examples

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

Derive the system:

.. code:: python

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
constants and specified quantities. Here, we specify sinusoidal forcing:

.. code:: python

   from numpy import array, linspace, sin
   from pydy.system import System

   sys = System(kane,
                constants={mass: 1.0, stiffness: 1.0,
                           damping: 0.2, gravity: 9.8},
                specifieds={force: lambda x, t: sin(t)},
                initial_conditions={position: 0.1, speed: -1.0},
                times=linspace(0.0, 10.0, 1000))

Integrate the equations of motion to get the state trajectories:

.. code:: python

   y = sys.integrate()

Plot the results:

.. code:: python

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

- nose_: 1.3.7
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
   (pydy-dev)$ pip install numpy scipy cython nose theano sympy ipython "notebook<5.0" "ipywidgets<5.0" version_information
   (pydy-dev)$ pip install matplotlib # make sure to do this after numpy
   (pydy-dev)$ git clone git@github.com:pydy/pydy.git
   (pydy-dev)$ cd pydy
   (pydy-dev)$ python setup.py develop

.. _virtualenvwrapper: https://pypi.python.org/pypi/virtualenvwrappe://pypi.python.org/pypi/virtualenvwrapper

Or with conda_::

   $ conda create -c pydy -n pydy-dev setuptools numpy scipy ipython "notebook<5.0" "ipywidgets<5.0" cython nose theano sympy matplotlib version_information
   $ source activate pydy-dev
   (pydy-dev)$ git clone git@github.com:pydy/pydy.git
   (pydy-dev)$ cd pydy
   (pydy-dev)$ conda develop .

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

0.5.0
-----

- SymPy introduced a backward incompatibility to differentiation Matrices in
  SymPy 1.2, which remained in SymPy 1.3, see:
  https://github.com/sympy/sympy/issues/14958. This breaks PyDy's System class,
  see: https://github.com/pydy/pydy/issues/395. A fix is introduced to handle
  all support versions of SymPy. [PR `#408`_]
- Added a new example for anthropomorphic arm. [PR `#406`_]
- Fixed errors in the differential drive example. [PR `#405`_]
- Added a new example for a scara arm. [PR `#402`_]
- Fixed errors due to backwards incompatible changes with various dependencies. [PR `#397`_]
- ODEFunctionGenerator now works with no constants symbols. [PR `#391`_]

.. _#408: https://github.com/pydy/pydy/pull/408
.. _#406: https://github.com/pydy/pydy/pull/406
.. _#405: https://github.com/pydy/pydy/pull/405
.. _#402: https://github.com/pydy/pydy/pull/402
.. _#397: https://github.com/pydy/pydy/pull/397
.. _#391: https://github.com/pydy/pydy/pull/391

0.4.0
-----

- Bumped minimum Jupyter notebook to 4.0 and restricted to < 5.0. [PR `#381`_]
- Removed several deprecated functions. [PR `#375`_]
- Bumped minimum required hard dependencies to Ubuntu 16.04 LTS package
  versions. [PR `#372`_]
- Implemented ThreeJS Tube Geometry. [PR `#368`_]
- Improved circle rendering. [PR `#357`_]
- kwargs can be passed from System.generate_ode_function to the matrix
  generator. [PR `#356`_]
- Lagrangian simple pendulum example added. [PR `#351`_]
- Derivatives can now be used as specifies in System. [PR `#340`_]
- The initial conditions can now be adjusted in the notebook GUI. [PR `#333`_]
- The width of the viz canvas is now properly bounded in the notebook. [PR `#332`_]
- Planes now render both sides in the visualization GUI. [PR `#330`_]
- Adds in more type checks for System.times. [PR `#322`_]
- Added an OctaveMatrixGenerator for basic Octave/Matlab printing. [PR `#323`_]
- Simplified the right hand side evaluation code in the ODEFunctionGenerator.
  Note that this change comes with some performance hits. [PR `#301`_]

.. _#381: https://github.com/pydy/pydy/pull/381
.. _#375: https://github.com/pydy/pydy/pull/375
.. _#372: https://github.com/pydy/pydy/pull/372
.. _#368: https://github.com/pydy/pydy/pull/368
.. _#357: https://github.com/pydy/pydy/pull/357
.. _#356: https://github.com/pydy/pydy/pull/356
.. _#351: https://github.com/pydy/pydy/pull/351
.. _#340: https://github.com/pydy/pydy/pull/340
.. _#333: https://github.com/pydy/pydy/pull/333
.. _#332: https://github.com/pydy/pydy/pull/332
.. _#330: https://github.com/pydy/pydy/pull/330
.. _#322: https://github.com/pydy/pydy/pull/322
.. _#323: https://github.com/pydy/pydy/pull/323
.. _#301: https://github.com/pydy/pydy/pull/301

0.3.1
-----

- Removed the general deprecation warning from System. [PR `#262`_]
- Don't assume user enters input in server shutdown. [PR `#264`_]
- Use vectorized operations to compute transformations. [PR `#266`_]
- Speedup theano generators. [PR `#267`_]
- Correct time is displayed on the animation slider. [PR `#272`_]
- Test optional dependencies only if installed. [PR `#276`_]
- Require benchmark to run in Travis. [PR `#277`_]
- Fix dependency minimum versions in setup.py [PR `#279`_]
- Make CSE optional in CMatrixGenerator. [PR `#284`_]
- Fix codegen line break. [PR `#292`_]
- Don't assume Scene always has a System. [PR `#295`_]
- Python 3.5 support and testing against Python 3.5 on Travis. [PR `#305`_]
- Set minimum dependency versions to match Ubuntu Trusty 14.04 LTS. [PR `#306`_]
- Replace sympy.phyics.mechanics deprecated methods. [PR `#309`_]
- Updated installation details to work with IPython/Jupyter 4.0. [PR `#311`_]
- Avoid the IPython widget deprecation warning if possible. [PR `#311`_]
- Updated the mass-spring-damper example to IPy4 and added version_information. [PR `#312`_]
- The Cython backend now compiles on Windows. [PR `#313`_]
- CI testing is now run on appveyor with Windows VMs. [PR `#315`_]
- Added a verbose option to the Cython compilation. [PR `#315`_]
- Fixed the RHS autogeneration. [PR `#318`_]
- Improved the camera code through inheritance [PR `#319`_]

.. _#262: https://github.com/pydy/pydy/pull/262
.. _#264: https://github.com/pydy/pydy/pull/264
.. _#266: https://github.com/pydy/pydy/pull/266
.. _#267: https://github.com/pydy/pydy/pull/267
.. _#272: https://github.com/pydy/pydy/pull/272
.. _#276: https://github.com/pydy/pydy/pull/276
.. _#277: https://github.com/pydy/pydy/pull/277
.. _#279: https://github.com/pydy/pydy/pull/279
.. _#284: https://github.com/pydy/pydy/pull/284
.. _#292: https://github.com/pydy/pydy/pull/292
.. _#295: https://github.com/pydy/pydy/pull/295
.. _#305: https://github.com/pydy/pydy/pull/305
.. _#306: https://github.com/pydy/pydy/pull/306
.. _#309: https://github.com/pydy/pydy/pull/309
.. _#311: https://github.com/pydy/pydy/pull/311
.. _#312: https://github.com/pydy/pydy/pull/312
.. _#313: https://github.com/pydy/pydy/pull/313
.. _#315: https://github.com/pydy/pydy/pull/315
.. _#318: https://github.com/pydy/pydy/pull/318
.. _#319: https://github.com/pydy/pydy/pull/319

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
