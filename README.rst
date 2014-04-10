====
PyDy
====

.. image:: https://travis-ci.org/pydy/pydy.png?branch=master
   :target: https://travis-ci.org/pydy/pydy

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
   >>> from pydy import codegen, viz

Viz
===
pydy-viz
========

Visualization of multibody systems generated with PyDy.

Installation
============

This package relies on dependencies included in the SciPy stack (i.e. NumPy and
matplotlib). These packages are not necessarily easy to install from source, so
it is best to follow the instructions available on the `SciPy installation
page`_.

.. _SciPy installation page: http://www.scipy.org/install.html

One example of installing the setuptools, NumPy, and matplotlib dependencies
for Debian based Linux systems is to install from the apt package manager::

   $ apt-get install python-setuptools python-numpy python-matplotlib

Once the dependencies are installed, then download the source and install with::

   $ python setup.py install

This will automatically install the latest version of the final dependency
SymPy if needed.

**Note that pydy-viz is currently only tested on Python 2.7.**

Tests
=====

The Python tests require nose so get them with your package manager::

   $ apt-get python-nose python-coverage

or pip::

   $ pip install nose coverage

The tests can be run from the root directory with::

   $ nosetests

And to see more detail with coverage, run::

   $ nosetests -v --with-coverage --cover-package=pydy_viz

These are alternative ways to run the Python tests::

   $ bin/test
   $ python setup.py nosetests

For the Javascript tests the Jasmine and blanket.js libraries are used.  Both
of these libraries are included in pydy-viz with the source. To run the
Javascript tests, go to the javascript library directory::

   $ cd pydy_viz/static/js

Then run a simple HTTP Server with Python (the server is required due to some
cross browser issues with blanket.js)::

   $ python -m SimpleHTTPServer

Now visit http://localhost:8000/SpecRunner.html in a webgl compliant browser.

Documentation
=============

The documentation is hosted at http://pydy-viz.readthedocs.org but you can also
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

Release Notes
=============

0.1.0
-----

- Initial release.

Code Generation
===============

This package provides code generation facilities for PyDy_. For now, it
generates functions that can evaluate the right hand side of the ordinary
differential equations generated with sympy.physics.mechanics_ with three
different backends: SymPy's lambdify_, Theano_, and Cython_.

.. _PyDy: http://pydy.org
.. _sympy.physics.mechanics: http://docs.sympy.org/latest/modules/physics/mechanics
.. _lambdify: http://docs.sympy.org/latest/modules/utilities/lambdify.html#sympy.utilities.lambdify.lambdify
.. _Theano: http://deeplearning.net/software/theano/
.. _Cython: http://cython.org/

Dependencies
------------

Required
~~~~~~~~

- Python: 2.7 (Python 3+ may work)
- setuptools
- NumPy: >=1.6.1
- SymPy: >=0.7.3

Optional
~~~~~~~~

- Cython: >=0.15.1
- Theano: >=0.6.0
- SciPy: >=0.9 (only for full examples)
- matplotlib: >=0.99 (only for full examples)

There are a variety of methods to install these packages. Refer to the SciPy
Stack installation instructions for details.

Installation
------------

Once the dependencies are installed, the package can be installed from PyPi
using::

   $ easy_install pydy-code-gen

or::

   $ pip install pydy-code-gen

You can also grab the source and then install\ [#]_.

Using the zip download::

   $ wget https://github.com/PythonDynamics/pydy-code-gen/archive/master.zip
   $ unzip pydy-code-gen-master.zip
   $ cd pydy-code-gen-master
   $ python setup.py install

Using Git::

   $ git clone https://github.com/PythonDynamics/pydy-code-gen.git
   $ cd pydy-code-gen
   $ python setup.py install

.. [#] Note that this is the latest development version. Specific releases
   can be found here: https://github.com/PythonDynamics/pydy-code-gen/releases
   or by checking out a tag with Git.

Usage
-----

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

Store the expressions and symbols in sequences for the code generation::

   mass_matrix = kane.mass_matrix_full
   forcing_vector = kane.forcing_full
   constants = (mass, stiffness, damping, gravity)
   coordinates = (position,)
   speeds = (speed,)
   specified = (force,)

Now generate the function needed for numerical evaluation of the ODEs. The
generator can use various back ends: ``lambdify``, ``theano``, or ``cython``::

   from pydy_code_gen.code import generate_ode_function

   evaluate_ode = generate_ode_function(mass_matrix, forcing_vector, constants,
                                        coordinates, speeds, specified,
                                        generator='lambdify')

Integrate the equations of motion under the influence of a specified sinusoidal
force::

   from numpy import array, linspace, sin
   from scipy.integrate import odeint

   x0 = array([0.1, -1.0])
   args = {'constants': array([1.0, 1.0, 0.2, 9.8]),
           'specified': lambda x, t: sin(t)}
   t = linspace(0.0, 10.0, 1000)

   y = odeint(evaluate_ode, x0, t, args=(args,))

Plot the results::

   import matplotlib.pyplot as plt

   plt.plot(t, y)
   plt.legend((str(position), str(speed)))
   plt.show()

Development Environment
-----------------------

Development Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

- nose: 1.3.0

Installation
~~~~~~~~~~~~

The following installation assumes you have virtualenvwrapper_ and all the
dependencies needed to build the packages::

   $ mkvirtualenv pydy-dev
   (pydy-dev)$ pip install numpy scipy cython nose theano sympy
   (pydy-dev)$ pip install matplotlib # make sure to do this after numpy
   (pydy-dev)$ git clone git@github.com:PythonDynamics/pydy-code-gen.git
   (pydy-dev)$ cd pydy-code-gen
   (pydy-dev)$ python setup.py develop

.. _virtualenvwrapper: https://pypi.python.org/pypi/virtualenvwrappe://pypi.python.org/pypi/virtualenvwrapper

Run the tests::

   (pydy-dev)$ nosetests

Run the benchmark to test the n-link pendulum problem.::

   (pydy-dev)$ python bin/benchmark_pydy_code_gen.py <max # of links> <# of time steps>

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
