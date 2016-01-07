#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages

exec(open('pydy/version.py').read())

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
    install_requires=['numpy>=1.8.1',
                      'scipy>=0.13.3',
                      'sympy>=0.7.4.1',
                      ],
    extras_require={'doc': ['sphinx', 'numpydoc'],
                    'codegen': ['Cython>=0.20.1',
                                'Theano>=0.7.0'],
                    'examples': ['matplotlib>=0.99',
                                 'ipython[notebook]>=3.0.0'],
                    },
    tests_require=['nose>=1.3.0'],
    test_suite='nose.collector',
    include_package_data=True,
    classifiers=['Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
)
