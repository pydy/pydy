PyDy Code Generation
====================

This distribution provides code generation facilities for PyDy_. For now, it
generates functions that can evaluate the right hand side of the ordinary
differential equations generated with sympy.physics.mechanics_ with three
different backends: SymPy's lambdify_, Theano_, and Cython_.

.. _PyDy: http://pydy.org
.. _sympy.physics.mechanics: http://docs.sympy.org/latest/modules/physics/mechanics
.. _lambdify: http://docs.sympy.org/latest/modules/utilities/lambdify.html#sympy.utilities.lambdify.lambdify
.. _Theano: http://deeplearning.net/software/theano/
.. _Cython: http://cython.org/

Dependencies
============

Required
--------

- Python: 2.7 (Python 3+ may work)
- setuptools
- NumPy: >=1.6.1
- SymPy: >=0.7.3

Optional
--------

- Cython: >=0.15.1
- Theano: >=0.6.0
- SciPy: >=0.9 (only for full examples)
- matplotlib: >=0.99 (only for full examples)

There are a variety of methods to install these packages. Refer to the SciPy
Stack installation instructions for details.

Installation
============

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
=======================

Development Dependencies
------------------------

- nose: 1.3.0

Installation
------------

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
