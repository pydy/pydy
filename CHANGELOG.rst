=============
Release Notes
=============

0.8.0 (TBA)
===========

0.7.1 (March 4, 2023)
=====================

- Reduced sdist size by moving the MANIFEST.in prune command last.

0.7.0 (March 4, 2023)
=====================

- Support Python 3.10 and 3.11. [PR `#488`_]
- Fixed the Carvallo-Whipple bicycle model to match Basu-Mandal benchmark
  numbers. [PR `#486`_]
- Added Box geometry to the javascript GUI [PR `#484`_]
- Updated the three link conical pendulum example to use new kanes_equations()
  syntax. [PR `#481`_]
- Added example of a 3D multilink pendulum with colliding bobs. [PR `#467`_]
- ``LambdifyODEFunctionGenerator`` now accepts a ``cse=True/False`` kwarg and
  if SymPy >=1.9 is installed, then the underlying generated code by
  ``lambdify`` will be simplified. It is ``True`` by default. [PR `#464`_]
- Visualization supports durations that don't start at 0.

.. _#464: https://github.com/pydy/pydy/pull/464
.. _#467: https://github.com/pydy/pydy/pull/467
.. _#481: https://github.com/pydy/pydy/pull/481
.. _#484: https://github.com/pydy/pydy/pull/484
.. _#486: https://github.com/pydy/pydy/pull/486
.. _#488: https://github.com/pydy/pydy/pull/488

0.6.0 (February 4, 2022)
========================

- Dropped support for Python 2.7 and 3.6. [PR `#459`_]
- Moved chaos pendulum example to Sphinx docs.
- Added Astrobee example [PR `#453`_]
- Added the ability to pass optional arguments to the ODE solver in System. [PR
  `#447`_]
- Cylinders, Spheres, and Circles loaded via PyThreeJS will appear more round.
  [PR `#440`_]
- Added a Carvallo-Whipple bicycle example to the documentation [PR `#442`_]
- Oldest supported dependencies for Python 3 are aligned with Ubuntu 20.04 LTS.
  For Python 2, the oldest necessary dependencies are used if the ones for
  Ubuntu 20.04 LTS are too new. [PR `#432`_]
- Dropped support for Python 3.5 [PR `#429`_]
- Improved the README and documentation integration. [PR `#424`_]
- Moved some examples to Sphinx [PRs `#421`_, `#423`_]
- jupyter-sphinx enabled for examples in the documentation [PR `#419`_]
- Added an example with no constraints that uses ``display_jupyter()`` for
  animation. [PR `#418`_]
- Added an example that has both configuration and motion constraints.
  [PR `#417`_]
- ``display_jupyter()`` method added to ``Scene`` that utilizes pythreejs for
  animating a system. [PR `#416`_]
- Remove support for required dependencies prior to those in Ubuntu 18.04 LTS.
  [PR `#415`_]
- Recommend installing from Conda Forge [PR `#411`_]

.. _#459: https://github.com/pydy/pydy/pull/459
.. _#453: https://github.com/pydy/pydy/pull/453
.. _#447: https://github.com/pydy/pydy/pull/447
.. _#442: https://github.com/pydy/pydy/pull/442
.. _#440: https://github.com/pydy/pydy/pull/440
.. _#432: https://github.com/pydy/pydy/pull/432
.. _#429: https://github.com/pydy/pydy/pull/429
.. _#424: https://github.com/pydy/pydy/pull/424
.. _#423: https://github.com/pydy/pydy/pull/423
.. _#421: https://github.com/pydy/pydy/pull/421
.. _#419: https://github.com/pydy/pydy/pull/419
.. _#418: https://github.com/pydy/pydy/pull/418
.. _#417: https://github.com/pydy/pydy/pull/417
.. _#416: https://github.com/pydy/pydy/pull/416
.. _#415: https://github.com/pydy/pydy/pull/415
.. _#411: https://github.com/pydy/pydy/pull/411

0.5.0 (January 9, 2019)
=======================

- SymPy introduced a backward incompatibility to differentiation Matrices in
  SymPy 1.2, which remained in SymPy 1.3, see:
  https://github.com/sympy/sympy/issues/14958. This breaks PyDy's System class,
  see: https://github.com/pydy/pydy/issues/395. A fix is introduced to handle
  all support versions of SymPy. [PR `#408`_]
- Added a new example for anthropomorphic arm. [PR `#406`_]
- Fixed errors in the differential drive example. [PR `#405`_]
- Added a new example for a scara arm. [PR `#402`_]
- Fixed errors due to backwards incompatible changes with various dependencies. [PR `#397`_]
- ODEFunctionGenerator now works with no constants symbols. [PR `#391`_]

.. _#408: https://github.com/pydy/pydy/pull/408
.. _#406: https://github.com/pydy/pydy/pull/406
.. _#405: https://github.com/pydy/pydy/pull/405
.. _#402: https://github.com/pydy/pydy/pull/402
.. _#397: https://github.com/pydy/pydy/pull/397
.. _#391: https://github.com/pydy/pydy/pull/391

0.4.0 (May 30, 2017)
====================

- Bumped minimum Jupyter notebook to 4.0 and restricted to < 5.0. [PR `#381`_]
- Removed several deprecated functions. [PR `#375`_]
- Bumped minimum required hard dependencies to Ubuntu 16.04 LTS package
  versions. [PR `#372`_]
- Implemented ThreeJS Tube Geometry. [PR `#368`_]
- Improved circle rendering. [PR `#357`_]
- kwargs can be passed from System.generate_ode_function to the matrix
  generator. [PR `#356`_]
- Lagrangian simple pendulum example added. [PR `#351`_]
- Derivatives can now be used as specifies in System. [PR `#340`_]
- The initial conditions can now be adjusted in the notebook GUI. [PR `#333`_]
- The width of the viz canvas is now properly bounded in the notebook. [PR `#332`_]
- Planes now render both sides in the visualization GUI. [PR `#330`_]
- Adds in more type checks for System.times. [PR `#322`_]
- Added an OctaveMatrixGenerator for basic Octave/Matlab printing. [PR `#323`_]
- Simplified the right hand side evaluation code in the ODEFunctionGenerator.
  Note that this change comes with some performance hits. [PR `#301`_]

.. _#381: https://github.com/pydy/pydy/pull/381
.. _#375: https://github.com/pydy/pydy/pull/375
.. _#372: https://github.com/pydy/pydy/pull/372
.. _#368: https://github.com/pydy/pydy/pull/368
.. _#357: https://github.com/pydy/pydy/pull/357
.. _#356: https://github.com/pydy/pydy/pull/356
.. _#351: https://github.com/pydy/pydy/pull/351
.. _#340: https://github.com/pydy/pydy/pull/340
.. _#333: https://github.com/pydy/pydy/pull/333
.. _#332: https://github.com/pydy/pydy/pull/332
.. _#330: https://github.com/pydy/pydy/pull/330
.. _#322: https://github.com/pydy/pydy/pull/322
.. _#323: https://github.com/pydy/pydy/pull/323
.. _#301: https://github.com/pydy/pydy/pull/301

0.3.1 (January 6, 2016)
=======================

- Removed the general deprecation warning from System. [PR `#262`_]
- Don't assume user enters input in server shutdown. [PR `#264`_]
- Use vectorized operations to compute transformations. [PR `#266`_]
- Speedup theano generators. [PR `#267`_]
- Correct time is displayed on the animation slider. [PR `#272`_]
- Test optional dependencies only if installed. [PR `#276`_]
- Require benchmark to run in Travis. [PR `#277`_]
- Fix dependency minimum versions in setup.py [PR `#279`_]
- Make CSE optional in CMatrixGenerator. [PR `#284`_]
- Fix codegen line break. [PR `#292`_]
- Don't assume Scene always has a System. [PR `#295`_]
- Python 3.5 support and testing against Python 3.5 on Travis. [PR `#305`_]
- Set minimum dependency versions to match Ubuntu Trusty 14.04 LTS. [PR `#306`_]
- Replace sympy.phyics.mechanics deprecated methods. [PR `#309`_]
- Updated installation details to work with IPython/Jupyter 4.0. [PR `#311`_]
- Avoid the IPython widget deprecation warning if possible. [PR `#311`_]
- Updated the mass-spring-damper example to IPy4 and added version_information. [PR `#312`_]
- The Cython backend now compiles on Windows. [PR `#313`_]
- CI testing is now run on appveyor with Windows VMs. [PR `#315`_]
- Added a verbose option to the Cython compilation. [PR `#315`_]
- Fixed the RHS autogeneration. [PR `#318`_]
- Improved the camera code through inheritance [PR `#319`_]

.. _#262: https://github.com/pydy/pydy/pull/262
.. _#264: https://github.com/pydy/pydy/pull/264
.. _#266: https://github.com/pydy/pydy/pull/266
.. _#267: https://github.com/pydy/pydy/pull/267
.. _#272: https://github.com/pydy/pydy/pull/272
.. _#276: https://github.com/pydy/pydy/pull/276
.. _#277: https://github.com/pydy/pydy/pull/277
.. _#279: https://github.com/pydy/pydy/pull/279
.. _#284: https://github.com/pydy/pydy/pull/284
.. _#292: https://github.com/pydy/pydy/pull/292
.. _#295: https://github.com/pydy/pydy/pull/295
.. _#305: https://github.com/pydy/pydy/pull/305
.. _#306: https://github.com/pydy/pydy/pull/306
.. _#309: https://github.com/pydy/pydy/pull/309
.. _#311: https://github.com/pydy/pydy/pull/311
.. _#312: https://github.com/pydy/pydy/pull/312
.. _#313: https://github.com/pydy/pydy/pull/313
.. _#315: https://github.com/pydy/pydy/pull/315
.. _#318: https://github.com/pydy/pydy/pull/318
.. _#319: https://github.com/pydy/pydy/pull/319

0.3.0 (January 19, 2015)
========================

User Facing
-----------

- Introduced conda builds and binstar support. [PR `#219`_]
- Dropped support for IPython < 3.0. [PR `#237`_]
- Added support Python 3.3 and 3.4. [PR `#229`_]
- Bumped up the minimum dependencies for NumPy, SciPy, and Cython [PR `#233`_].
- Removed the partial implementation of the Mesh shape. [PR `#172`_]
- Overhauled the code generation package to make the generators more easily
  extensible and to improve simulation speed. [PR `#113`_]
- The visualizer has been overhauled as part of Tarun Gaba's 2014 GSoC
  internship [PR `#82`_]. Here are some of the changes:

  - The JavaScript is now handled by AJAX and requires a simple server.
  - The JavaScript has been overhauled and now uses prototype.js for object
    oriented design.
  - The visualizer can now be loaded in an IPython notebook via IPython's
    widgets using ``Scene.display_ipython()``.
  - A slider was added to manually control the frame playback.
  - The visualization shapes' attributes can be manipulated via the GUI.
  - The scene json file can be edited and downloaded from the GUI.
  - pydy.viz generates two JSONs now (instead of one in earlier versions). The
    JSON generated from earlier versions will **not** work in the new version.
  - Shapes can now have a material attribute.
  - Model constants can be modified and the simulations can be rerun all via
    the GUI.
  - Switched from socket based server to python's core SimpleHTTPServer.
  - The server has a proper shutdown response [PR `#241`_]

- Added a new experimental System class and module to more seamlessly manage
  integrating the equations of motion. [PR `#81`_]

.. _#241: https://github.com/pydy/pydy/pull/241
.. _#237: https://github.com/pydy/pydy/pull/237
.. _#229: https://github.com/pydy/pydy/pull/229
.. _#233: https://github.com/pydy/pydy/pull/233
.. _#219: https://github.com/pydy/pydy/pull/219
.. _#172: https://github.com/pydy/pydy/pull/172
.. _#113: https://github.com/pydy/pydy/pull/113
.. _#82: https://github.com/pydy/pydy/pull/82
.. _#81: https://github.com/pydy/pydy/pull/81

Development
-----------

- Switched to a conda based Travis testing setup. [PR `#231`_]
- When using older SymPy development versions with non-PEP440 compliant version
  identifiers, setuptools < 8 is required. [PR `#166`_]
- Development version numbers are now PEP 440 compliant. [PR `#141`_]
- Introduced pull request checklists and CONTRIBUTING file. [PR `#146`_]
- Introduced light code linting into Travis. [PR `#148`_]

.. _#231: https://github.com/pydy/pydy/pull/231
.. _#166: https://github.com/pydy/pydy/pull/166
.. _#141: https://github.com/pydy/pydy/pull/141
.. _#146: https://github.com/pydy/pydy/pull/146
.. _#148: https://github.com/pydy/pydy/pull/148

0.2.1 (June 19, 2014)
=====================

- Unbundled unnecessary files from tar ball.

0.2.0 (June 19, 2014)
=====================

- Merged pydy_viz, pydy_code_gen, and pydy_examples into the source tree.
- Added a method to output "static" visualizations from a Scene object.
- Dropped the matplotlib dependency and now only three.js colors are valid.
- Added joint torques to the n_pendulum model.
- Added basic examples for codegen and viz.
- Graceful fail if theano or cython are not present.
- Shapes can now use sympy symbols for geometric dimensions.
