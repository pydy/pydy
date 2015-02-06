#!/usr/bin/env python

from ..utils import sympy_equal_to_or_newer_than


def test_sympy_equal_to_or_newer_than():

    assert sympy_equal_to_or_newer_than('0.7.6-git', '0.7.6-git')
    assert sympy_equal_to_or_newer_than('0.7.6', '0.7.6-git')
    assert sympy_equal_to_or_newer_than('0.7.5', '0.7.6-git')
    assert sympy_equal_to_or_newer_than('0.6.5', '0.7.6-git')
    assert not sympy_equal_to_or_newer_than('0.7.7', '0.7.6-git')
