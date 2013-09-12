Installation
------------


There are different methods you can use to get PyDyViz up and running for your system.

Dependencies:

1) Python >=2.5
2) SymPy >=0.7.2
3) NumPy >=
4) SciPy >=

SciPy is used by the Code Generator module for numerical integration of Equations of Motions.
So not exactly a dependency for PyDyViz, but is required anyways for Code Generator to work.

Source
======

PyDyViz can be installed from source via archive. Download the latest release archive from here.
Extract the archive and cd into the directory. 
Issue the following command from command line::

    $ python setup.py install

It should get the PyDyViz up and running.
        

Git
===

If you are a developer, or looking for the latest features without waiting for another release, you can install the development version 
from the git repository.

Issue the following command from terminal::

    $ git clone git://github.com/PythonDynamics/pydy-viz.git

and then change directory to the cloned directory and run setup.py install::

    $ cd pydy-viz
    $ python setup.py install
    
and you will have the latest development version on board.
    


Python Package Index
====================

You can also install PyDyViz from Python Package Index using pip or easy_install

Its as easy as opening a terminal and typing in::

    $ pip install pydy_viz
    
or for easy_install::

    $ easy_install pydy_viz   

    
Run PyDyViz
===========

You can check the installation by running the following command from Python Interpreter::

    >>> import pydy_viz

If it does not throws any error/traceback, it means that PyDyViz is installed.    


You can check the version of the software by issuing following command from interpreter::

    >>> import pydy_viz
    >>> print pydy_viz.__version__

Questions
=========

If you have any question about installation, or any general question, feel free to visit the IRC channel at irc.freenode.net, channel #pydy. 
In addition, our mailing list is an excellent source of community support.

If you think there’s a bug or you would like to request a feature, please open an issue ticket.

