#!/usr/bin/env python

from pkg_resources import parse_version
from setuptools import __version__ as SETUPTOOLS_VERSION
from nose.tools import assert_raises

from ..utils import sympy_equal_to_or_newer_than


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
