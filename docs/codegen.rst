=======
codegen
=======

Introduction
============

The :py:mod:`pydy.codegen` package contains various tools to generate numerical
code from symbolic descriptions of the equations of motion of systems. It
allows you to generate code using a variety of backends depending on your
needs. The generated code can also be auto-wrapped for immediate use in a
Python session or script. Each component of the code generators and wrappers
are accessible so that you can use just the raw code or the wrapper versions.

We currently support three backends:

``lambdify``
   This generates NumPy-aware Python code which is defined in a Python
   ``lambda`` function, using the ``sympy.utilities.lambdify`` module and is
   the default generator.
``Theano``
   This generates Theano trees that are compiled into low level code, using the
   ``sympy.printers.theano_code`` module.
``Cython``
   This generates C code that can be called from Python, using SymPy's C code
   printer utilities and Cython.

On Windows
==========

For the Cython backend to work on Windows you must install a suitable compiler.
See this `Cython wiki page
<https://github.com/cython/cython/wiki/CythonExtensionsOnWindows>`_ for
instructions on getting a compiler installed. The easiest solution is to use
the Microsoft Visual C++ Compiler for Python.

Example Use
===========

The simplest entry point to the code generation tools is through the
:py:class:`~pydy.system.System` class.

.. code:: pycon

   >>> from pydy.models import multi_mass_spring_damper
   >>> sys = multi_mass_spring_damper()
   >>> type(sys)
   <class 'pydy.system.System'>
   >>> rhs = sys.generate_ode_function()
   >>> help(rhs) # rhs is a function:
   Returns the derivatives of the states, i.e. numerically evaluates the right
   hand side of the first order differential equation.

   x' = f(x, t, p)

   Parameters
   ==========
   x : ndarray, shape(2,)
       The state vector is ordered as such:
           - x0(t)
           - v0(t)
   t : float
       The current time.
   p : dictionary len(3) or ndarray shape(3,)
       Either a dictionary that maps the constants symbols to their numerical
       values or an array with the constants in the following order:
           - m0
           - c0
           - k0

   Returns
   =======
   dx : ndarray, shape(2,)
       The derivative of the state vector.


   >>> import numpy as np
   >>> rhs(np.array([1.0, 2.0]), 0.0, np.array([1.0, 2.0, 3.0]))
   array([ 2., -7.])

You can also use the functional interface to the code generation/wrapper
classes:

.. code:: pycon

   >>> from numpy import array
   >>> from pydy.models import multi_mass_spring_damper
   >>> from pydy.codegen.ode_function_generators import generate_ode_function
   >>> sys = multi_mass_spring_damper()
   >>> sym_rhs = sys.eom_method.rhs()
   >>> q = sys.coordinates
   >>> u = sys.speeds
   >>> p = sys.constants_symbols
   >>> rhs = generate_ode_function(sym_rhs, q, u, p)
   >>> rhs(array([1.0, 2.0]), 0.0, array([1.0, 2.0, 3.0]))
   array([ 2., -7.])

Other backends can be used by passing in the ``generator`` keyword argument,
e.g.:

.. code:: pycon

   >>> rhs = generate_ode_function(sym_rhs, q, u, p, generator='cython')
   >>> rhs(array([1.0, 2.0]), 0.0, array([1.0, 2.0, 3.0]))
   array([ 2., -7.])

The backends are implemented as subclasses of
:py:class:`~pydy.codegen.ode_function_generators.ODEFunctionGenerator`. You can
make use of the ``ODEFunctionGenerator`` classes directly:

.. code:: pycon

   >>> from pydy.codegen.ode_function_generators import LambdifyODEFunctionGenerator
   >>> g = LambdifyODEFunctionGenerator(sym_rhs, q, u, p)
   >>> rhs = g.generate()
   >>> rhs(array([1.0, 2.0]), 0.0, array([1.0, 2.0, 3.0]))
   array([ 2., -7.])

Furthermore, for direct control over evaluating matrices you can use the
``lamdify`` and ``theano_functions`` in SymPy or utilize the
:py:class:`~pydy.codegen.cython_code.CythonMatrixGenerator` class
in PyDy. For example, this shows you how to generate C and Cython code to
evaluate matrices:

.. code:: pycon

   >>> from pydy.codegen.cython_code import CythonMatrixGenerator
   >>> sys = multi_mass_spring_damper()
   >>> q = sys.coordinates
   >>> u = sys.speeds
   >>> p = sys.constants_symbols
   >>> sym_rhs = sys.eom_method.rhs()
   >>> g = CythonMatrixGenerator([q, u, p], [sym_rhs])
   >>> setup_py, cython_src, c_header, c_src = g.doprint()
   >>> print(setup_py)
   #!/usr/bin/env python

   from distutils.core import setup
   from distutils.extension import Extension

   from Cython.Build import cythonize
   import numpy

   extension = Extension(name="pydy_codegen",
                         sources=["pydy_codegen.pyx",
                                  "pydy_codegen_c.c"],
                         include_dirs=[numpy.get_include()])

   setup(name="pydy_codegen",
         ext_modules=cythonize([extension]))

   >>> print(cython_src)
   import numpy as np
   cimport numpy as np
   cimport cython

   cdef extern from "pydy_codegen_c.h":
       void evaluate(
                     double* input_0,
                     double* input_1,
                     double* input_2,
                     double* output_0
                    )

   @cython.boundscheck(False)
   @cython.wraparound(False)
   def eval(
            np.ndarray[np.double_t, ndim=1, mode='c'] input_0,
            np.ndarray[np.double_t, ndim=1, mode='c'] input_1,
            np.ndarray[np.double_t, ndim=1, mode='c'] input_2,
            np.ndarray[np.double_t, ndim=1, mode='c'] output_0
           ):

       evaluate(
                <double*> input_0.data,
                <double*> input_1.data,
                <double*> input_2.data,
                <double*> output_0.data
               )

       return (
               output_0
              )

   >>> print(c_src)
   #include <math.h>
   #include "pydy_codegen_c.h"

   void evaluate(
                 double input_0[1],
                 double input_1[1],
                 double input_2[3],
                 double output_0[2]
                )
   {

       double pydy_0 = input_1[0];

       output_0[0] = pydy_0;
       output_0[1] = (-input_2[1]*pydy_0 - input_2[2]*input_0[0])/input_2[0];

   }

   >>> print(c_header)
   void evaluate(
                 double input_0[1],
                 double input_1[1],
                 double input_2[3],
                 double output_0[2]
                );
   /*

   input_0[1] : [x0(t)]
   input_1[1] : [v0(t)]
   input_2[3] : [m0, c0, k0]

   */

   >>> rhs = g.compile()
   >>> res = array([0.0, 0.0])
   >>> rhs(array([1.0]), array([2.0]), array([1.0, 2.0, 3.0]), res)
   array([ 2., -7.])

We also support generating Octave/Matlab code as shown below:

.. code:: pycon

   >>> from pydy.codegen.octave_code import OctaveMatrixGenerator
   >>> sys = multi_mass_spring_damper()
   >>> q = sys.coordinates
   >>> u = sys.speeds
   >>> p = sys.constants_symbols
   >>> sym_rhs = sys.eom_method.rhs()
   >>> g = OctaveMatrixGenerator([q + u, p], [sym_rhs])
   >>> m_src = g.doprint()
   >>> print(m_src)
   function [output_1] = eval_mats(input_1, input_2)
   % function [output_1] = eval_mats(input_1, input_2)
   %
   % input_1 : [x0(t), v0(t)]
   % input_2 : [k0, m0, c0]

       pydy_0 = input_1(2);

       output_1 = [pydy_0; (-input_2(3).*pydy_0 - ...
       input_2(1).*input_1(1))./input_2(2)];

   end

API
===

.. automodule:: pydy.codegen.c_code
   :members:
   :special-members: __init__

.. automodule:: pydy.codegen.cython_code
   :members:
   :special-members: __init__

.. automodule:: pydy.codegen.matrix_generator
   :members:
   :special-members: __init__

.. automodule:: pydy.codegen.octave_code
   :members:
   :special-members: __init__

.. automodule:: pydy.codegen.ode_function_generators
   :members:
   :special-members: __init__
