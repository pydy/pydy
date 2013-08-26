#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='PyDy',
    version='0.1.0dev',
    author='PyDy Authors',
    author_email='pydy@googlegroups.com',
    url="http://pydy.org",
    description='Python tool kit for multi-body dynamics.',
    long_description=open('README.rst').read(),
    keywords="multibody dynamics",
    license='LICENSE.txt',
    packages=find_packages(),
    install_requires=['sympy>=0.7.2',
                      'pydy-viz>=0.1.0',
                      'numpy',
                      'scipy',
                      'ipython[notebook]',
                      ],
    include_package_data=True,
    classifiers=['Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
)
