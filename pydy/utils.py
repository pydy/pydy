#!/usr/bin/env python

import textwrap

from pkg_resources import parse_version
import sympy as sm

SYMPY_VERSION = sm.__version__


def sympy_equal_to_or_newer_than(version):
    """Returns true if the installed version of SymPy is equal to or newer
    than the provided version string."""
    return cmp(parse_version(SYMPY_VERSION), parse_version(version)) > -1


def wrap_and_indent(lines, indentation=4, width=79):
    """Returns a single string in which the lines have been indented and
    wrapped into a block of text."""
    # TODO : This will indent any lines that only contain a new line. Which
    # may not be preferable.
    new_lines = []
    for line in lines:
        if line != '\n':
            wrapped = textwrap.wrap(line, width=width-indentation)
        else:
            wrapped = [line]
        new_lines += wrapped
    spacer = '\n' + ' ' * indentation
    return ' ' * indentation + spacer.join(new_lines)
