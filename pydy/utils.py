#!/usr/bin/env python

from pkg_resources import parse_version
import sympy as sm

SYMPY_VERSION = sm.__version__


def sympy_equal_to_or_newer_than(version):
    """Returns true if the installed version of SymPy is equal to or newer
    than the provided version string."""
    return cmp(parse_version(SYMPY_VERSION), parse_version(version)) > -1
