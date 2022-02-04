#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages

exec(open('pydy/version.py').read())

# The lower bounds for the dependencies match those in Ubuntu 20.04 LTS.

install_requires = [
    'numpy>=1.16.5',
    'setuptools>=44.0.0',
    'scipy>=1.3.3',
    'sympy>=1.5.1',
]

extras_require = {
    'doc': [
        'jupyter_sphinx',  # not in Ubuntu repos
        'numpydoc>=0.70',
        'pythreejs',  # not in Ubuntu repos
        'sphinx>=1.8.5',
    ],
    'codegen': [
        'Cython>=0.29.14',
        'Theano>=1.0.4'
    ],
    'examples': [
        'matplotlib>=3.1.2',
        'notebook>=4.0.0,<5.0.0',  # for display_ipython()
        'ipywidgets>=4.0.0,<5.0.0'  # for display_ipython()
    ],
}

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
    extras_require=extras_require,
    tests_require=['nose>=1.3.7'],
    test_suite='nose.collector',
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
