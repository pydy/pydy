#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = '0.1.0dev'

try:
    import pydy_code_gen.code as codegen
except ImportError:
    print('pydy-code-gen not installed.')

try:
    import pydy_viz as visualization
except ImportError:
    print('pydy-viz not installed.')
