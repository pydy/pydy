====
PyDy
====

.. image:: https://travis-ci.org/pydy/pydy.png?branch=master
   :target: https://travis-ci.org/pydy/pydy

PyDy_, short for Python Dynamics, is a tool kit written in and accessed through
the Python programming language that utilizes an array of scientific tools to
study multibody dynamics. The goal is to have a modular framework and
eventually a physics abstraction layer which utilizes a variety of backends
that can provide the user with their desired workflow, including:

- Model specification
- Equation of motion generation
- Simulation
- Visualization
- Publication

We started by building the SymPy_ `mechanics package`_ which provides an API
for building models and generating the symbolic equations of motion for complex
multibody systems and have more recently developed two packages, pydy.codegen
and pydy.viz, for simulation and visualization of the models. The remaining
tools currently used in the PyDy workflow are popular scientific Python
packages such as NumPy_, SciPy_, IPython_, and matplotlib_ (i.e. the SciPy
stack) which provide additional code for numerical analyses, simulation, and
visualization.

Installation
============

The PyDy workflow has hard dependencies on these Python packages:

- Python >= 2.7
- setuptools

SciPy Stack

- SymPy_ >= 0.7.4.1
- NumPy_ >= 1.6.1
- SciPy_ >= 0.9.0
- IPython_ (3.0.0-dev)

It's best to install the SciPy Stack dependencies using the instructions_
provided on the SciPy website.

**Note:** The pydy.viz has patched IPython(3.0.0-dev) as the dependency. It can be installed by::

    $ pip install -e git+https://github.com/ipython/ipython.git


Once the dependencies are installed, the package can be installed from PyPi
using::

   $ easy_install pydy

or::

   $ pip install pydy

For system wide installs you will need root permissions (perhaps prepend
commands with ``sudo``).

You can also grab the source and then install\ [#]_.

Using the zip download::

   $ wget https://github.com/pydy/pydy/archive/master.zip
   $ unzip pydy-master.zip
   $ cd pydy-master
   $ python setup.py install

Using Git::

   $ git clone https://github.com/pydy/pydy.git
   $ cd pydy
   $ python setup.py install

.. [#] Note that this is the latest development version. Specific releases
   can be found here: https://github.com/pydy/pydy/releases
   or by checking out a tag with Git.

Development Environment
-----------------------

Development Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

Tests require nose:

- nose: 1.3.0

Javascript testing requires phantomjs:

- phantomjs: 1.9.0

Isolated Environments
~~~~~~~~~~~~~~~~~~~~~

The following installation assumes you have virtualenvwrapper_ and all the
dependencies needed to build the various packages::

   $ mkvirtualenv pydy-dev
   (pydy-dev)$ pip install numpy scipy cython nose theano sympy
   (pydy-dev)$ pip install matplotlib # make sure to do this after numpy
   (pydy-dev)$ pip install -e git+https://github.com/ipython/ipython.git
   (pydy-dev)$ git clone git@github.com:pydy/pydy.git
   (pydy-dev)$ cd pydy
   (pydy-dev)$ python setup.py develop

.. _virtualenvwrapper: https://pypi.python.org/pypi/virtualenvwrappe://pypi.python.org/pypi/virtualenvwrapper

Or with conda_::

   $ conda create -n pydy-dev numpy scipy cython nose theano sympy matplotlib
   $ source activate pydy-dev
   (pydy-dev)$ conda remove ipython  # conda automatically installs ipython into envs
   (pydy-dev)$ git clone https://github.com/ipython/ipython
   (pydy-dev)$ cd ipython
   (pydy-dev)$ git submodule update
   (pydy-dev)$ pip install -e ".[notebook]"
   (pydy-dev)$ cd ..
   (pydy-dev)$ git clone git@github.com:pydy/pydy.git
   (pydy-dev)$ cd pydy
   (pydy-dev)$ python setup.py develop

.. _conda: https://github.com/conda/conda

Run the tests::

   (pydy-dev)$ nosetests

For the Javascript tests the Jasmine and blanket.js libraries are used. Both
of these libraries are included in pydy.viz with the source. To run the
Javascript tests::

   cd pydy/viz/static/js/tests && phantomjs run-jasmine.js SpecRunner.html && cd ../../../../../


Run the benchmark to test the n-link pendulum problem.::

   (pydy-dev)$ python bin/benchmark_pydy_code_gen.py <max # of links> <# of time steps>

Usage
=====

This is an example of a simple 1 degree of freedom system: a mass, spring,
damper system under the influence of gravity and a force::


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

Create a system to manage integration. Specify numerical values for the
constants and specified quantities. Here, we specify sinusoidal forcing::

   from numpy import array, linspace, sin
   from pydy.system import System

   sys = System(kane,
                constants={mass: 1.0, stiffness: 1.0,
                           damping: 0.2, gravity: 9.8},
                specified={'symbols': [force],
                           'values': lambda x, t: sin(t)},
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
build them from source using the following instructions:

Requires:

- Sphinx
- numpydoc

::

   pip install sphinx numpydoc

To build the HTML docs::

   $ sphinx-build -b html docs/src docs/build

View::

   $ firefox docs/build/index.html

Packages
========

Code Generation
---------------

This package provides code generation facilities for PyDy_. For now, it
generates functions that can evaluate the right hand side of the ordinary
differential equations generated with sympy.physics.mechanics_ with three
different backends: SymPy's lambdify_, Theano_, and Cython_.

.. _PyDy: http://pydy.org
.. _sympy.physics.mechanics: http://docs.sympy.org/latest/modules/physics/mechanics
.. _lambdify: http://docs.sympy.org/latest/modules/utilities/lambdify.html#sympy.utilities.lambdify.lambdify
.. _Theano: http://deeplearning.net/software/theano/
.. _Cython: http://cython.org/

To enable different code generation backends, you can install the various
optional dependencies:

- Cython: >=0.15.1
- Theano: >=0.6.0

Visualization (viz)
-------------------

Visualization of multibody systems generated with PyDy.

Related Packages
================

These are various related Python packages that have similar functionality.

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

Release Notes
=============

0.3.0
-----

- Added a new System class and module to more seamlessly manage integrating the
  equations of motion.
- Updated UI for pydy.viz's Visualizer.
- pydy.viz generates two JSONs now (instead of one in earlier versions). The JSON generated from earlier versions
  will **not** work in the new version.
- Upgraded IPython version dependency to '3.0.0-dev'. 
- Added new method `Scene.display_ipython` to allow visualizations to run in IPython notebooks (requires IPython 3.0.0-dev for this to work)
- Added new material attribute for Shape.
- Added IPython widgets, to change simulation parameters from the GUI itself.
- Now EoMs can be re-integrated by changing simulation parameters from GUI (by hitting `rerun simulation` button).
- Switched from socket based server to python's core SimpleHTTPServer.
- Automated javascript testing on travis.


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

Citation
========

If you make use of the PyDy toolchain in your work or research, please cite us
in your publications or on the web. This citation can be used:

   Gilbert Gede, Dale L Peterson, Angadh S Nanjangud, Jason K Moore, and Mont
   Hubbard, "Constrained Multibody Dynamics With Python: From Symbolic Equation
   Generation to Publication", ASME 2013 International Design Engineering
   Technical Conferences and Computers and Information in Engineering
   Conference, 2013, `10.1115/DETC2013-13470
   <http://dx.doi.org/10.1115/DETC2013-13470>`_.
