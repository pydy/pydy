#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages

exec(open('pydy/version.py').read())

# The lower bounds for the dependencies match those in Ubuntu 24.04 LTS.

install_requires = [
    'numpy>=1.26.4',
    'packaging>=24.0',
    'scipy>=1.11.4',
    'setuptools>=68.1.2',
    'sympy>=1.12',
]

extras_require = {
    'doc': [
        'jupyter_sphinx>=0.3.2',
        'numpydoc>=1.6.0',
        'pythreejs>=2.4.2',  # not in Ubuntu repos, 2023-02-20 PyPi release
        'sphinx>=7.2.6',
    ],
    'codegen': [
        'Cython>=0.29.37',  # cython3-legacy
        'Theano>=1.0.5',  # not in Ubuntu repos, 2020-07-27 PyPi release
        'symjit>=2.5.0',  # not in Ubuntu repos
    ],
    'examples': [
        'ipywidgets>=8.1.1',  # for display_ipython()
        'matplotlib>=3.6.3',
        'notebook>=6.4.12',  # for display_ipython()
    ],
}

if os.name == 'nt':
    install_requires.append('PyWin32>=306')  # 2023-03-26 PyPi release

setup(
    name='pydy',
    version=__version__,
    author='PyDy Authors',
    author_email='pydy@googlegroups.com',
    url="https://pydy.org",
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
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
