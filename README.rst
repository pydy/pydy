pydy-viz
========

Exploration of visualization for PyDy systems.

Installation
============

Probably easiest to install the main dependencies from your package manager,
e.g.:

$ apt-get python-numpy python-matplotlib

Download the source and then install with setuptools (this will pull in the
latest version of SymPy).

$ python setup.py install

Tests
=====

The Python tests require nose:


$ pip install nose coverage

nosetests -v --with-coverage --cover-package=pydy_viz

or run

bin/test

after nose is installed.
