#!/usr/bin/env python

import os

import numpy as np
import sympy as sm

from ...models import multi_mass_spring_damper
from ..cython_code import CythonMatrixGenerator


class TestCythonMatrixGenerator(object):

    def setup(self):

        self.prefix = 'boogly_bee'

        sys = multi_mass_spring_damper(6, True, True)

        self.matrices = (sys.eom_method.mass_matrix,
                         sys.eom_method.forcing)

        self.arguments = (sys.constants_symbols,
                          sys.coordinates,
                          sys.speeds,
                          sys.specifieds_symbols)

        self.generator = CythonMatrixGenerator(self.arguments,
                                               self.matrices, self.prefix)

    def test_generate_code_blocks(self):

        expected = {}

        expected['header_args'] = \
"""\
                  double* input_0,
                  double* input_1,
                  double* input_2,
                  double* input_3,
                  double* output_0,
                  double* output_1\
"""

        expected['python_args'] = \
"""\
         np.ndarray[np.double_t, ndim=1, mode='c'] input_0,
         np.ndarray[np.double_t, ndim=1, mode='c'] input_1,
         np.ndarray[np.double_t, ndim=1, mode='c'] input_2,
         np.ndarray[np.double_t, ndim=1, mode='c'] input_3,
         np.ndarray[np.double_t, ndim=1, mode='c'] output_0,
         np.ndarray[np.double_t, ndim=1, mode='c'] output_1\
"""

        expected['c_args'] = \
"""\
             <double*> input_0.data,
             <double*> input_1.data,
             <double*> input_2.data,
             <double*> input_3.data,
             <double*> output_0.data,
             <double*> output_1.data\
"""

        expected['output'] = \
"""\
            output_0.reshape(6, 6),
            output_1\
"""

        self.generator._generate_code_blocks()

        for k, v in self.generator.code_blocks.items():
            assert v == expected[k]

    def test_doprint(self):

        expected_pyx_source = \
"""\
import numpy as np
cimport numpy as np
cimport cython

cdef extern from "boogly_bee_c.h":
    void evaluate(
                  double* input_0,
                  double* input_1,
                  double* input_2,
                  double* input_3,
                  double* output_0,
                  double* output_1
                 )

@cython.boundscheck(False)
@cython.wraparound(False)
def eval(
         np.ndarray[np.double_t, ndim=1, mode='c'] input_0,
         np.ndarray[np.double_t, ndim=1, mode='c'] input_1,
         np.ndarray[np.double_t, ndim=1, mode='c'] input_2,
         np.ndarray[np.double_t, ndim=1, mode='c'] input_3,
         np.ndarray[np.double_t, ndim=1, mode='c'] output_0,
         np.ndarray[np.double_t, ndim=1, mode='c'] output_1
        ):

    evaluate(
             <double*> input_0.data,
             <double*> input_1.data,
             <double*> input_2.data,
             <double*> input_3.data,
             <double*> output_0.data,
             <double*> output_1.data
            )

    return (
            output_0.reshape(6, 6),
            output_1
           )\
"""

        expected_setup_py_source = """\
#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize
import numpy

extension = Extension(name="boogly_bee",
                      sources=["boogly_bee.pyx",
                               "boogly_bee_c.c"],
                      include_dirs=[numpy.get_include()])

setup(name="boogly_bee",
      ext_modules=cythonize([extension]))\
"""

        setup, pyx, c_header, c_source = self.generator.doprint()

        assert setup == expected_setup_py_source
        assert pyx == expected_pyx_source

    def test_write(self):

        setup, pyx, c_header, c_source = self.generator.doprint()

        self.generator.write()

        with open(self.prefix + '_c.h') as f:
            assert f.read() == c_header

        with open(self.prefix + '_c.c') as f:
            assert f.read() == c_source

        with open(self.prefix + '_setup.py') as f:
            assert f.read() == setup

        with open(self.prefix + '.pyx') as f:
            assert f.read() == pyx

    def test_compile(self):

        f = self.generator.compile()

        subs = {}

        args = []
        for argset in self.arguments:
            vals = np.random.random(len(argset))
            args.append(vals)
            for arg, val in zip(argset, vals):
                subs[arg] = val

        for matrix in self.matrices:
            nr, nc = matrix.shape
            args.append(np.empty(nr * nc, dtype=float))

        for output, expected in zip(f(*args), self.matrices):
            try:
                expected = sm.matrix2numpy(expected.subs(subs),
                                           dtype=float).squeeze()
            except TypeError:
                # dtype kwarg in not supported in earlier SymPy versions
                expected = np.asarray(sm.matrix2numpy(expected.subs(subs)),
                                      dtype=float).squeeze()

            np.testing.assert_allclose(output, expected)

    def teardown(self):

        for suffix in ['_c.h', '_c.c', '_setup.py', '.pyx']:
            filename = self.prefix + suffix
            if os.path.isfile(filename):
                os.remove(filename)
