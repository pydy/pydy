======
models
======

Introduction
============

The :py:mod:`pydy.models` file provides canned symbolic models of classical
dynamic systems that are mostly for testing and example purposes. There are
currently two models:

:py:func:`~pydy.models.multi_mass_spring_damper`
   A one dimensional series of masses connected by linear dampers and springs
   that can optionally be under the influence of gravity and an arbitrary
   force.
:py:func:`~pydy.models.n_link_pendulum_on_cart`
   This is an extension to the classic two dimensional inverted pendulum on a
   cart to multiple links. You can optionally apply an arbitrary lateral
   force to the cart and/or apply arbitrary torques between each link.

Example Use
===========

A simple one degree of freedom mass spring damper system can be created with:

.. code:: python

   >>> from pydy.models import multi_mass_spring_damper
   >>> sys = multi_mass_spring_damper()
   >>> sys.constants_symbols
   {m0, c0, k0}
   >>> sys.coordinates
   [x0(t)]
   >>> sys.speeds
   [v0(t)]
   >>> sys.eom_method.rhs()
   Matrix([
   [                    v0(t)],
   [(-c0*v0(t) - k0*x0(t))/m0]])

A two degree of freedom mass spring damper system under the influence of
gravity and two external forces can be created with:

.. code:: python

   >>> sys = multi_mass_spring_damper(2, True, True)
   >>> sys.constants_symbols
   {c1, m1, k0, c0, k1, m0, g}
   >>> sys.coordinates
   [x0(t), x1(t)]
   >>> sys.speeds
   [v0(t), v1(t)]
   >>> sys.specifieds_symbols
   {f0(t), f1(t)}
   >>> from sympy import simplify
   >>> sm.simplify(sys.eom_method.rhs())
   Matrix([
   [                                                                                                              v0(t)],
   [                                                                                                              v1(t)],
   [                                                     (-c0*v0(t) + c1*v1(t) + g*m0 - k0*x0(t) + k1*x1(t) + f0(t))/m0],
   [-(m1*(-c0*v0(t) + g*m0 + g*m1 - k0*x0(t) + f0(t) + f1(t)) + (m0 + m1)*(c1*v1(t) - g*m1 + k1*x1(t) - f1(t)))/(m0*m1)]])

API
===

.. automodule:: pydy.models
   :members:
   :special-members: __init__
