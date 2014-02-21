#!/usr/bin/env python

from setuptools import setup, find_packages

import pydy_viz

setup(
    name='pydy-viz',
    version=pydy_viz.__version__,
    author='Tarun Gaba',
    author_email='tarun.gaba7@gmail.com',
    url="https://github.com/pydy/pydy-viz/",
    description='Browser based 3D visualization of multibody simulations.',
    long_description=open('README.rst').read(),
    keywords="dynamics multibody simulation visualization",
    license='LICENSE.txt',
    packages=find_packages(),
    install_requires=['sympy>=0.7.2',
                      'numpy>=1.6.1',
                      'matplotlib>=0.99.0'],
    extras_require={'doc': ['sphinx', 'numpydoc']},
    tests_require=['nose'],
    test_suite='nose.collector',
    scripts=['bin/test'],
    include_package_data=True,
    classifiers=['Development Status :: 3 - Alpha',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
)
