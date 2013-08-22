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

The Python tests require nose and the javascript tests require Jasmine and
PhantomJS.

$ pip install nose coverage

nosetests -v --with-coverage --cover-package=pydy_viz

$ sudo add-apt-repository ppa:chris-lea/node.js
$ sudo apt-get update
$ sudo apt-get install nodejs npm phantomjs
$ sudo npm install -g phantom-jasmine
$ # may have to create this symlink
$ sudo ln -s /usr/bin/nodejs /usr/bin/node
$ phantom-jasmine pydy_viz/static/js/SpecRunner.html

or run

bin/test

after all is installed for both of those commands.
