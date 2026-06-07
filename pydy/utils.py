#!/usr/bin/env python

import re
import textwrap
import itertools

from packaging.version import parse as parse_version
import sympy as sm
import sympy.physics.mechanics as me
import numpy as np

SYMPY_VERSION = sm.__version__


def _sort_velocity_constraints(velocity_constraints, q, qdot_exprs):
    """Return the time differentiated holonomic constraints and nonholonomic
    constraints separately.

    Given velcoity leval constraings that may be a combination of time
    differentiated holonomic constraints and nonholonomic constraints, this
    function uses the symmetry of second derivatives check to distinguish the
    constraints.

    If the velocity constraints are large expressions, this could be slow due
    to the need to simplify for zero checking.

    Parameters
    ==========
    velocity_constraints : iterable of Expr
    q : iterable of Function of time
        Generalized coordinates.
    qdot_exprs : iteraable of Expr
        Expressions for the q' definitions.

    Returns
    =======
    holonomic_idxs : list of integer
        Indices that identify the rows of ``velocity_constraints`` which are
        time differentiated holonomic constraints.
    nonholonomic_idxs : list of integer
        Indices that identify the rows of ``velocity_constraints`` which are
        nonholonomic constraints.

    """
    # TODO : This could be something useful to implement in SymPy.
    velocity_constraints = sm.Matrix(velocity_constraints)
    jac, _ = sm.linear_eq_to_matrix(velocity_constraints, qdot_exprs)
    nonholonomic_idxs = []
    idxs = list(range(len(qdot_exprs)))
    comb_idxs = itertools.combinations(idxs, 2)
    for i, row in enumerate(jac.tolist()):
        for (j, k) in comb_idxs:
            zero = row[j].diff(q[k]) - row[k].diff(q[j])
            syms = list(zero.atoms(sm.Symbol) | me.find_dynamicsymbols(zero))
            eval_zero = sm.lambdify(syms, zero)
            vals = np.random.random(len(syms))
            if not np.allclose(eval_zero(*vals), 0.0, atol=1e-12):
                nonholonomic_idxs.append(i)
                break
    all_idxs = range(len(velocity_constraints))
    holonomic_idxs = [i for i in all_idxs if i not in nonholonomic_idxs]

    return holonomic_idxs, nonholonomic_idxs


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
