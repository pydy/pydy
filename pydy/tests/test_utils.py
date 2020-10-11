#!/usr/bin/env python

from pkg_resources import parse_version
from setuptools import __version__ as SETUPTOOLS_VERSION
from nose.tools import assert_raises

from sympy import cos, sin, tan, sqrt, Matrix
from sympy.physics.mechanics import dynamicsymbols

from ..utils import sympy_equal_to_or_newer_than, wrap_and_indent
from ..codegen.cython_code import CythonMatrixGenerator


def test_sympy_equal_to_or_newer_than():

    # sympy_equal_to_or_newer_than(version, installed_version)
    assert sympy_equal_to_or_newer_than('0.7.6.dev', '0.7.6.dev')
    assert not sympy_equal_to_or_newer_than('0.7.6', '0.7.6.dev')
    assert sympy_equal_to_or_newer_than('0.7.5', '0.7.6.dev')
    assert sympy_equal_to_or_newer_than('0.6.5', '0.7.6.dev')
    assert not sympy_equal_to_or_newer_than('0.7.7', '0.7.6.dev')
    if parse_version(SETUPTOOLS_VERSION) >= parse_version('8.0'):
        with assert_raises(ValueError):
            sympy_equal_to_or_newer_than('0.7.7', '0.7.6-git')


def test_codegen_linewrap():

    # Generated can result in long expressions with no obvious place to insert a
    # newline. Refer to https://github.com/pydy/pydy/issues/263.
    x, y, z = dynamicsymbols('x y z')
    expr = (x*y*z*cos(x)*sin(x)*tan(x)*sqrt(x)*cos(y)*sin(y)*tan(y)*
            sqrt(y)*cos(z)*sin(z)*tan(z)*sqrt(z)*cos(x*y)*sin(x*y)*tan(x*y)*
            sqrt(x*y)* cos(y*z)*sin(y*z)*tan(y*z)*sqrt(y*z)*cos(z*x)*
            sin(z*x)*tan(z*x)*sqrt(z*x)*3)/8
    mat_expr = Matrix([expr])

    q = [x, y, z]
    # Don't raise an error with this line.
    gen = CythonMatrixGenerator([q], [mat_expr])


def test_wrap_and_indent():

    lines = ["a + b + c + d + e", "a + b + c + d + e"]
    wrapped = wrap_and_indent(lines, width=10)
    expected = """\
    a + b
    + c +
    d + e
    a + b
    + c +
    d + e"""
    assert wrapped == expected

    wrapped = wrap_and_indent(lines, width=14, continuation=' ...')
    expected = """\
    a + b ...
    + c + ...
    d + e
    a + b ...
    + c + ...
    d + e"""
    assert wrapped == expected

    lines = ["% a + b + c + d + e", "% a + b + c + d + e"]
    wrapped = wrap_and_indent(lines, width=12, comment='%')
    expected = """\
    % a + b
    % + c + d
    % + e
    % a + b
    % + c + d
    % + e"""
    assert wrapped == expected
