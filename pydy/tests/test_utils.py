#!/usr/bin/env python

from sympy import Symbol

from pkg_resources import parse_version
from setuptools import __version__ as SETUPTOOLS_VERSION
from nose.tools import assert_raises

from ..bodies import Body
from ..utils import sympy_equal_to_or_newer_than, Mass, Length, \
        GravitationalConstant, AssignSymbol


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

def test_mass():

    mass = Mass()
    assert mass.symbol == Symbol('m')
    assert mass.value == 1

    body = Body('body')
    body_mass = Mass(body.name)
    assert body_mass.symbol == Symbol(body.name + "_mass")

    body_mass_value = Mass(body.name, 4)
    assert body_mass_value.value == 4

def test_spring_constant():

    spring_constant = SpringConstant()
    assert spring_constant.constant == Symbol('k')
    assert spring_constant.value == 1

    body = Body('body')
    body_mass = Mass(body.name)
    assert body_mass.symbol == Symbol(body.name + "_mass")

    body_mass_value = Mass(body.name, 4)
    assert body_mass_value.value == 4

def test_length():

    length = Length()
    assert length.constant == Symbol('l')
    assert length.value == 1

    body = Body('body')
    body_length = Length(body.name)
    assert body_length.symbol == Symbol(body.name + "_length")

    body_length_value = Length(body.name, 2)
    assert body_length_value.value == 2

def test_gravitational_constant():

    gravitational_constant = GravitationalConstant()
    assert gravitational_constant.constant == Symbol('g')
    assert gravitational_constant.value == 9.8

    body_gravity_value = GravitationalConstant(10)
    assert body_gravity_value.value == 10

def test_assign_symbol():
    assign_symbol = AssignSymbol()

    assert assign_symbol.get_body_name() == 'body_1'
    assert assign_symbol.get_joint_name() == 'joint_1'
    assert assign_symbol.get_length_name() == 'length_1'
    assert assign_symbol.get_mass_name() == 'mass_1'
    assert assign_symbol.get_symbol('coordinate') == 'coordinate_1'
    assert assign_symbol.get_symbol('coordinate') == 'coordinate_2'
