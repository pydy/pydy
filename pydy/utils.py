#!/usr/bin/env python

import re
import textwrap

from pkg_resources import parse_version
from setuptools import __version__ as SETUPTOOLS_VERSION
import sympy as sm
from sympy.core.function import AppliedUndef
from sympy.utilities.iterables import iterable
from sympy.physics.mechanics import dynamicsymbols

SYMPY_VERSION = sm.__version__


def sympy_equal_to_or_newer_than(version, installed_version=None):
    """Returns true if the installed version of SymPy is equal to or newer
    than the provided version string."""

    if installed_version is None:
        v = SYMPY_VERSION
    else:
        v = installed_version

    if v.endswith('-git') and \
            parse_version(SETUPTOOLS_VERSION) >= parse_version('8.0'):

        msg = ('You are using an older development version of SymPy with a '
               'non-PEP440 compliant version number: {}. Please install '
               'setuptools < 8.0 or a newer development version of SymPy.')
        raise ValueError(msg.format(v))

    return parse_version(v) >= parse_version(version)


def sympy_newer_than(version):
    """Returns true if the installed version of SymPy is newer than the
    provided version string."""
    return parse_version(SYMPY_VERSION) > parse_version(version)


def wrap_and_indent(lines, indentation=4, width=79, continuation=None,
                    comment=None):
    """Returns a single string in which the lines have been indented and
    wrapped into a block of text.

    Parameters
    ==========
    indentation : integer
        The number of characters to indent.
    width : integer
        The maximum line width.
    continuation : string
        The continuation characters.
    comment : string
        The character that designates a comment line.

    """

    if continuation is None:
        cont_len = 0
    else:
        cont_len = len(continuation)

    if comment is None:
        comment_len = 0
    else:
        comment_len = len(comment)

    # TODO : This will indent any lines that only contain a new line. Which
    # may not be preferable.
    new_lines = []

    # TODO : The Octave printer has ".*" and "./" as operators and this doesn't
    # deal with that.
    # add whitespace before and after [*/] binary operands between
    # subexpressions and input/output
    pattern = re.compile('(\w\])([*/])(\w)')
    for line in lines:
        if line != '\n':
            line = pattern.sub(lambda m: ' '.join(m.groups()), line)
            wrapped = textwrap.wrap(line,
                                    width=width-indentation-cont_len-comment_len,
                                    break_long_words=False)
            if continuation:
                last = wrapped[-1]
                wrapped = [l + continuation for l in wrapped[:-1]]
                wrapped.append(last)

            if comment:
                for i, l in enumerate(wrapped[1:]):
                    wrapped[i + 1] = comment + ' ' + l
        else:
            wrapped = [line]
        new_lines += wrapped
    spacer = '\n' + ' ' * indentation
    return ' ' * indentation + spacer.join(new_lines)


# This is a copy of the function in SymPy. It doesn't exist in SymPy 0.7.4
# so we keep it here for now.
def find_dynamicsymbols(expression, exclude=None):
    """Find all dynamicsymbols in expression.

    >>> from sympy.physics.mechanics import dynamicsymbols, find_dynamicsymbols
    >>> x, y = dynamicsymbols('x, y')
    >>> expr = x + x.diff()*y
    >>> find_dynamicsymbols(expr)
    set([x(t), y(t), Derivative(x(t), t)])

    If the optional ``exclude`` kwarg is used, only dynamicsymbols
    not in the iterable ``exclude`` are returned.

    >>> find_dynamicsymbols(expr, [x, y])
    set([Derivative(x(t), t)])
    """
    t_set = set([dynamicsymbols._t])
    if exclude:
        if iterable(exclude):
            exclude_set = set(exclude)
        else:
            raise TypeError("exclude kwarg must be iterable")
    else:
        exclude_set = set()
    return set([i for i in expression.atoms(AppliedUndef, sm.Derivative) if
                i.free_symbols == t_set]) - exclude_set


class PyDyDeprecationWarning(DeprecationWarning):
    pass

class PyDyImportWarning(ImportWarning):
    pass

class PyDyFutureWarning(FutureWarning):
    pass

class PyDyUserWarning(UserWarning):
    pass
