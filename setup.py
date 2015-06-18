#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

from setuptools import setup, find_packages

exec(open('pydy/version.py').read())

# I was getting the same error as:
# https://github.com/statsmodels/statsmodels/issues/1073, so the following
# line is added.
os.environ["MPLCONFIGDIR"] = "."

if sys.version_info >= (3, 0):
    NUMPY_MIN_VER = '1.9'
    SCIPY_MIN_VER = '0.14.0'
    SYMPY_MIN_VER = '0.7.5'
    CYTHON_MIN_VER = '0.20.1'
    THEANO_MIN_VER = '0.7.0'
else:
    NUMPY_MIN_VER = '1.7'
    SCIPY_MIN_VER = '0.11'
    SYMPY_MIN_VER = '0.7.4.1'
    CYTHON_MIN_VER = '0.17'
    THEANO_MIN_VER = '0.6.0'

setup(
    name='pydy',
    version=__version__,
    author='PyDy Authors',
    author_email='pydy@googlegroups.com',
    url="http://pydy.org",
    description='Python tool kit for multi-body dynamics.',
    long_description=open('README.rst').read(),
    keywords="multibody dynamics",
    license='LICENSE.txt',
    packages=find_packages(),
    install_requires=['numpy>={}'.format(NUMPY_MIN_VER),
                      'scipy>={}'.format(SCIPY_MIN_VER),
                      'sympy>={}'.format(SYMPY_MIN_VER),
                      ],
    extras_require={'doc': ['sphinx', 'numpydoc'],
                    'codegen': ['Cython>={}'.format(CYTHON_MIN_VER),
                                'Theano>={}'.format(THEANO_MIN_VER)],
                    'examples': ['matplotlib>=0.99',
                                 'ipython[notebook]>=0.3.0'],
                    },
    tests_require=['nose>=1.3.0'],
    test_suite='nose.collector',
    scripts=['bin/benchmark_pydy_code_gen.py'],
    include_package_data=True,
    classifiers=['Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
)
