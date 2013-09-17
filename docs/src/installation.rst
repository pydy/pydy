Installation
------------


There are different methods you can use to get PyDyViz up and running 
for your system.

Dependencies:

1) Python 
2) SymPy 
3) NumPy 

SciPy is used by the Code Generator module for numerical integration 
of Equations of Motions.So not exactly a dependency for PyDyViz, 
but is required anyways for Code Generator to work.

Git
===

If you are a developer, or looking for the latest features without 
waiting for another release, you can install the development version 
from the git repository.

Issue the following command from terminal::

    $ git clone git://github.com/PythonDynamics/pydy-viz.git

and then change directory to the cloned directory and run setup.py 
install::

    $ cd pydy-viz
    $ python setup.py install
    
and you will have the latest development version on board.
    
Run PyDyViz
===========

You can check the installation by running the following command from 
Python Interpreter::

    >>> import pydy_viz

If it does not throws any error/traceback, it means that PyDyViz 
is installed.    


You can check the version of the software by issuing following 
command from interpreter::

    >>> import pydy_viz
    >>> print pydy_viz.__version__

Questions
=========

If you have any question about installation, or any general question, 
feel free to visit the IRC channel at irc.freenode.net, channel `#pydy`_. 
In addition, our `mailing list`_ is an excellent source of 
community support.

If you think there’s a bug or you would like to request a feature, 
please open an `issue`_.


.. _issue: https://github.com/PythonDynamics/pydy-viz/issues
.. _mailing list: http://groups.google.com/group/pydy
.. _#pydy: irc://irc.freenode.net/pydy
