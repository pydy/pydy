====================
PyDy's Documentation
====================

.. include:: index-opening.txt

Table of Contents
=================

.. toctree::
   :maxdepth: 2

   install.rst
   usage.rst
   tutorials.rst
   codegen.rst
   models.rst
   system.rst
   viz/viz.rst
   viz/api.rst
   utils.rst
   changelog.rst

.. ifconfig:: INCLUDE_EXAMPLES

   .. toctree::
      :hidden:

      examples/astrobee.rst
      examples/carvallo-whipple.rst
      examples/chaos-pendulum.rst
      examples/kane-levinson-1985-chapter-02.rst
      examples/kane-levinson-1985-chapter-03.rst
      examples/mass-spring-damper.rst
      examples/multidof-holonomic.rst
      examples/rocket-car.rst
      examples/three-link-conical-pendulum.rst
      examples/3d-n-body-pendulum

.. ifconfig:: INCLUDE_EXAMPLES

   Examples
   ========

   .. list-table::

      * - .. figure:: examples/mass_spring_damper.svg
             :width: 200px
             :target: examples/mass-spring-damper.html

             Linear mass-spring-damper system with gravity.
        - .. figure:: examples/multidof-holonomic.png
             :width: 200px
             :target: examples/multidof-holonomic.html

             A double compound and simple pendulum.
      * - .. figure:: examples/three-link-conical-pendulum.gif
             :width: 200px
             :target: examples/three-link-conical-pendulum.html

             Three link conical compound pendulum.
        - .. figure:: examples/bicycle-geometry.png
             :width: 200px
             :target: examples/carvallo-whipple.html

             Carvallo-Whipple bicycle model.
      * - .. figure:: examples/astrobee.png
             :width: 200px
             :target: examples/astrobee.html

             Astrobee free-flying ISS robot.
        - .. figure:: examples/chaos-pendulum.svg
             :width: 200px
             :target: examples/chaos-pendulum.html

             3D perpendicular axis double pendulum that exhibits chaos.

      * - .. figure:: examples/kane-levinson-1985.png
             :width: 200px
             :target: examples/kane-levinson-1985-chapter-02.html

             Exercises from Chapter 2 of Kane & Levinson 1985.
        - .. figure:: examples/kane-levinson-1985.png
             :width: 200px
             :target: examples/kane-levinson-1985-chapter-03.html

             Exercises from Chapter 3 of Kane & Levinson 1985.
      * - .. figure:: examples/rocket-car.png
             :width: 200px
             :target: examples/rocket-car.html

             Nonholonomic Rocket Car.
        -

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
