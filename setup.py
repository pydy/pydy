#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages

exec(open('pydy/version.py').read())

# The lower bounds for the dependencies match those in Ubuntu 22.04 LTS.

install_requires = [
    'numpy>=1.21.5',
    'scipy>=1.8.0',
    'setuptools>=44.1.1',
    'sympy>=1.9',
]

extras_require = {
    'doc': [
        'jupyter_sphinx>=0.3.2',
        'numpydoc>=1.2',
        'pythreejs',  # not in Ubuntu repos
        'sphinx>=4.3.2',
    ],
    'codegen': [
        'Cython>=0.29.28',
        'Theano>=1.0.5'
    ],
    'examples': [
        'matplotlib>=3.5.1',
        'notebook>=4.0.0,<5.0.0',  # for display_ipython()
        'ipywidgets>=4.0.0,<5.0.0'  # for display_ipython()
    ],
}

if os.name == 'nt':
    install_requires.append('PyWin32>=303')

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
    tests_require=['pytest'],
    include_package_data=True,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
