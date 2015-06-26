#!/usr/bin/env python

SYMPY_VERSION = sm.__version__


class TestBodies():

    def test_body_init(self):
        body = Body('body')

        assert body.name is None
        assert body.mass == Symbol('mass')

        sys = System(method='joints')
        sys.add_body(body)

        assert body.name == 'body_1'

        body_length = Length(body.name)
        assert body_length.symbol == Symbol(body.name + "_length")

        body_length_value = Length(body.name, 2)
        assert body_length_value.value == 2
