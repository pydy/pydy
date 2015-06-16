#!/usr/bin/env python

"""This module remains for backwards compatibility reasons and will be
removed in PyDy 0.4.0."""

from pkg_resources import parse_version
import warnings

import pydy
from ..utils import PyDyDeprecationWarning
from .ode_function_generators import generate_ode_function as new_gen_ode_func


if parse_version(pydy.__version__) > parse_version('0.4.0'):
    msg = ("This module, 'pydy.codegen.code', is no longer supported. "
           "Please remove this module.")
    raise ValueError(msg)

warnings.warn("This module, 'pydy.codgen.code', is deprecated. The "
              "function 'generate_ode_function' can be found in the "
              "'pydy.codegen.ode_function_generator' module. "
              "'CythonGenerator' has been removed, use "
              "'pydy.codegen.cython_code.CythonMatrixGenerator' "
              "instead.",
              PyDyDeprecationWarning)


class CythonGenerator(object):
    def __init__(self, *args, **kwargs):
        warnings.warn("'CythonGenerator' has been removed, use "
                      "'pydy.codegen.cython_code.CythonMatrixGenerator' "
                      "instead.", PyDyDeprecationWarning)


def generate_ode_function(mass_matrix, forcing_vector, constants,
                          coordinates, speeds, specified=None,
                          generator='lambdify'):
    """Returns a numerical function which can evaluate the right hand side
    of the first order ordinary differential equations from a system
    described by:

    M(constants, coordinates) x' = F(constants, coordinates, speeds, specified)

    Parameters
    ----------
    mass_matrix : sympy.Matrix, shape(n,n)
        The symbolic mass matrix of the system. The rows should correspond
        to the coordinates and speeds.
    forcing_vector : sympy.Matrix, shape(n,1)
        The symbolic forcing vector of the system.
    constants : list of sympy.Symbol
        The constants in the equations of motion.
    coordinates : list of sympy.Function
        The generalized coordinates of the system.
    speeds : list of sympy.Function
        The generalized speeds of the system.
    specified : list of sympy.Function
        The specifed quantities of the system.
    generator : string, {'lambdify'|'theano'|'cython'}, optional
        The method used for generating the numeric right hand side.

    Returns
    -------
    evaluate_ode_function : function
        A function which evaluates the derivaties of the states.

    """
    warnings.warn("This function is deprecated and will be removed in PyDy "
                  "0.4.0. Use the the new 'generate_ode_function' in "
                  "'pydy.codegen.ode_function_generator'",
                  PyDyDeprecationWarning)

    return new_gen_ode_func(forcing_vector, coordinates, speeds, constants,
                            mass_matrix=mass_matrix, specifieds=specified,
                            generator=generator)
