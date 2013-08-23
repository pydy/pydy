from setuptools import setup, find_packages

setup(
    name='pydy-viz',
    version='0.1.0dev',
    author='Tarun Gaba',
    author_email='tarun.gaba7@gmail.com',
    url="https://github.com/PythonDynamics/pydy-viz/",
    description='Browser based 3D visualization of multibody simulations.',
    long_description=open('README.rst').read(),
    keywords="dynamics multibody simulation visualization",
    license='LICENSE.txt',
    packages=find_packages(),
    install_requires=['sympy>=0.7.2', 'numpy'],#, 'matplotlib'],
    # For some reason matplotlib doesn't install if numpy isn't already
    # installed, even if pip is bringing in all the deps simultaneously. I
    # think matplotlib must not be using the setuptools "install_requires"
    # option and this balks. We will not actually depend on matplotlib, but
    # we do for now.
    extras_require={'doc': ['sphinx', 'numpydoc']},
    tests_require=['nose'],
    test_suite='nose.collector',
    scripts=['bin/test'],
    include_package_data=True,
    classifiers=[
                 'Development Status :: 2 - Pre-Alpha',
                 'Intended Audience :: Science/Research',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python :: 2.7',
                 'Topic :: Scientific/Engineering :: Physics',
                ],
)
