#!/usr/bin/env python

from setuptools import setup, find_packages

from pydy_code_gen import __version__

setup(
    name='pydy-code-gen',
    version=__version__,
    author='Jason K. Moore',
    author_email='moorepants@gmail.com',
    url="https://github.com/PythonDynamics/pydy-code-gen/",
    description='Code generation for multibody dynamic systems.',
    long_description=open('README.rst').read(),
    keywords="dynamics multibody simulation code generation",
    license='UNLICENSE',
    packages=find_packages(),
    install_requires=['sympy>=0.7.3',
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
    classifiers=[
                 'Development Status :: 3 - Alpha',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Physics',
    ],
)
