pydy-viz
========

Visualization of multibody systems generated with PyDy.

Installation
============

Probably easiest to install the main dependencies from your package manager,
e.g.::

   $ apt-get python-numpy python-setuptools

Then download the source and install with setuptools (this will automatically
install the latest version of SymPy)::

   $ python setup.py install

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
