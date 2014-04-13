import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_module = Extension(name="desired_mass_forcing",
                       sources=["desired_mass_forcing.pyx",
                                "desired_mass_forcing_c.c"],
                       include_dirs=[numpy.get_include()])

setup(name="desired_mass_forcing",
      cmdclass = {'build_ext': build_ext},
      ext_modules = [ext_module])
