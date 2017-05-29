#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages

exec(open('pydy/version.py').read())

# I was getting the same error as:
# https://github.com/statsmodels/statsmodels/issues/1073, so the following
# line is added.
os.environ["MPLCONFIGDIR"] = "."


install_requires=['numpy>=1.11.0',
                  'scipy>=0.17.1',
                  'sympy>=0.7.6.1']

if os.name == 'nt':
    install_requires.append('PyWin32>=219')

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
    install_requires=install_requires,
    extras_require={'doc': ['sphinx', 'numpydoc'],
                    'codegen': ['Cython>=0.23.4',
                                'Theano>=0.8.0'],
                    'examples': ['matplotlib>=1.5.1',
                                 'notebook>=4.0.0,<5.0.0',
                                 'ipywidgets>=4.0.0,<5.0.0'],
                    },
    tests_require=['nose>=1.3.7'],
    test_suite='nose.collector',
    include_package_data=True,
    classifiers=['Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3.5',
                 'Programming Language :: Python :: 3.6',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
)
