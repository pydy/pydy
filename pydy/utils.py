#!/usr/bin/env python

import re
import textwrap

from packaging.version import parse as parse_version
import sympy as sm

SYMPY_VERSION = sm.__version__


def sympy_equal_to_or_newer_than(version, installed_version=None):
    """Returns true if the installed version of SymPy is equal to or newer
    than the provided version string."""

    if installed_version is None:
        v = SYMPY_VERSION
    else:
        v = installed_version

    if v.endswith('-git'):
        msg = ('You are using an older development version of SymPy with a '
               'non-PEP440 compliant version number: {}. Please install '
               'a newer development version of SymPy.')
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
    pattern = re.compile(r'(\w\])([*/])(\w)')
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


class PyDyDeprecationWarning(DeprecationWarning):
    pass


class PyDyImportWarning(ImportWarning):
    pass


class PyDyFutureWarning(FutureWarning):
    pass


class PyDyUserWarning(UserWarning):
    pass
