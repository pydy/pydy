=============================
PyDy Package's documentation!
=============================

This is the central page for all PyDy's Documentation.

.. include:: ../README.rst

system module
=============

.. toctree::
   :maxdepth: 2

   system.rst

models module
=============

.. toctree::
   :maxdepth: 1

   models.rst

codegen package
===============

.. toctree::
    :maxdepth: 2

    codegen/codegen.rst
    codegen/api.rst

viz package
===========

.. toctree::
    :maxdepth: 2

    viz/viz.rst
    viz/using_visualizer.rst
    viz/api.rst

Examples
========

.. toctree::
    :maxdepth: 2

    tutorials/beginners.rst
    tutorials/advanced.rst

.. ifconfig:: INCLUDE_EXAMPLES

   .. list-table::

      * - .. figure:: examples/mass_spring_damper.svg
             :width: 200px
             :target: examples/mass-spring-damper.html

             Linear mass-spring-damper system with gravity.
        - .. figure:: examples/multidof-holonomic.png
             :width: 200px
             :target: examples/multidof-holonomic.html

             A double compound and simple pendulum.
        - .. figure:: examples/three-link-conical-pendulum.gif
             :width: 200px
             :target: examples/three-link-conical-pendulum.html

             Three link conical compound pendulum.
        -
      * - .. figure:: examples/kane-levinson-1985.png
             :width: 200px
             :target: examples/kane-levinson-1985-chapter-02.html

             Exercises from Chapter 2 of Kane & Levinson 1985.
        - .. figure:: examples/kane-levinson-1985.png
             :width: 200px
             :target: examples/kane-levinson-1985-chapter-03.html

             Exercises from Chapter 3 of Kane & Levinson 1985.
        -
        -

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
