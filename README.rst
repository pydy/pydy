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
