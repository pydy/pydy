#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages

from pydy import __version__

# I was getting the same error as:
# https://github.com/statsmodels/statsmodels/issues/1073, so the following
# line is added.
os.environ["MPLCONFIGDIR"] = "."

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
    install_requires=['sympy>=0.7.4.1',
                      'numpy>=1.6.1',
                      ],
    extras_require={'doc': ['sphinx', 'numpydoc'],
                    'theano': ['Theano>=0.6.0'],
                    'cython': ['Cython>=0.15.1'],
                    'examples': ['scipy>=0.9', 'matplotlib>=0.99'],
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
