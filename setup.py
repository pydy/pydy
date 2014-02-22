#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

# TODO : Figure out how to import the version from pydy/__init__.py even
# though pydy-code-gen and pydy-viz are not installed.

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
    install_requires=['pydy-viz>=0.1.0',
                      'pydy-code-gen>=0.1.0',
                      ],
    include_package_data=True,
    classifiers=['Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
)
