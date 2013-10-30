PyDy Code Gen
=============

This distribution provides some code generation facilities for PyDy. For now,
it generates functions that can evaluate the right hand side of the ordinary
differential equations generated with sympy.physics.mechanics with three
different backends: SymPy's lambdify, Theano, and Cython.

Dependencies
============

NumPy: 1.7.1
Cython: 0.19.2
SymPy: master
Theano: master

Developement Environmet Installation Process
============================================

Development Dependencies
------------------------

NumPy: 1.7.1
SciPy: 0.13.0
Cython: 0.19.2
nose: 1.3.0
matplotlib: 1.3.1
SymPy: master (>0.7.3)
Theano: master (>0.6.0rc3)

The following installation assumes you have virtualenv wrapper and all the
dependencies needed to build the packages::

   $ mkvirtualenv pydy-dev
   (pydy-dev)$ pip install numpy scipy cython nose
   (pydy-dev)$ pip install matplotlib # make sure to do this after numpy
   (pydy-dev)$ git clone theano
   (pydy-dev)$ cd theano
   (pydy-dev)$ python setup.py install
   (pydy-dev)$ cd ..
   (pydy-dev)$ git clone sympy
   (pydy-dev)$ cd sympy
   (pydy-dev)$ python setup.py install
   (pydy-dev)$ cd ..
   (pydy-dev)$ git clone pydy-code-gen
   (pydy-dev)$ cd pydy-code-gen
   (pydy-dev)$ python setup.py install
   (pydy-dev)$ cd misc
   (pydy-dev)$ python benchmark.py 5 1000
