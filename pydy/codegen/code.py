#!/usr/bin/env python

"""This module remains for backwards compatibility reasons and will be
removed in PyDy 0.4.0."""

import warnings

from .ode_function_generator import generate_ode_function

warnings.warn("This module, 'pydy.codgen.code', is deprecated. The function "
              "'generate_ode_function' can be found in the "
              "'pydy.codegen.ode_function_generator' module. "
              "'CythonGenerator' has been removed, use "
              "'pydy.codegen.cython_code.CythonMatrixGenerator' instead.",
              DeprecationWarning)
