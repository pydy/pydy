#!/usr/bin/env python

import os
import sys
import shutil
import tempfile
import importlib
import subprocess
from collections import defaultdict

from .c_code import CMatrixGenerator
from ..utils import wrap_and_indent


class CythonMatrixGenerator(object):
    """This class generates the Cython code for evaluating a sequence of
    matrices. It can compile the code and return a Python function."""

    _pyx_template = \
"""\
import numpy as np
cimport numpy as np
cimport cython

cdef extern from "{prefix}_c.h":
    void evaluate(
{header_args}
                 )

@cython.boundscheck(False)
@cython.wraparound(False)
def eval(
{python_args}
        ):

    evaluate(
{c_args}
            )

    return (
{output}
           )\
"""

    _setup_py_template = """\
#!/usr/bin/env python

from setuptools import setup
from setuptools import Extension

from Cython.Build import cythonize
import numpy

extension = Extension(name="{prefix}",
                      sources=["{prefix}.pyx",
                               "{prefix}_c.c"],
                      include_dirs=[numpy.get_include()])

setup(name="{prefix}",
      ext_modules=cythonize([extension]))\
"""

    _module_counter = 0

    def __init__(self, arguments, matrices, prefix='pydy_codegen', cse=True):
        """

        Parameters
        ==========
        arguments : sequences of sequences of SymPy Symbol or Function.
            Each of the sequences will be converted to input arrays in the
            Cython function. All of the symbols/functions contained in
            ``matrices`` need to be in the sequences, but the sequences can
            also contain extra symbols/functions that are not contained in
            the matrices.
        matrices : sequence of SymPy.Matrix
            A sequence of the matrices that should be evaluated in the
            function. The expressions should contain only sympy.Symbol or
            sympy.Function that are functions of me.dynamicsymbols._t.
        prefix : string, optional
            The desired prefix for the generated files.
        cse : boolean
            Find and replace common sub-expressions in ``matrices`` if True.

        """

        self.prefix = prefix
        self.matrices = matrices
        self.arguments = arguments
        self.num_matrices = len(matrices)
        self.num_arguments = len(arguments)
        self.c_matrix_generator = CMatrixGenerator(arguments, matrices,
                                                   cse=cse)

        self._generate_code_blocks()

    def _generate_code_blocks(self):

        lines = defaultdict(list)

        hd = 'double* {}_{},'
        py = "np.ndarray[np.double_t, ndim=1, mode='c'] {}_{},"
        c = '<double*> {}_{}.data,'
        out = 'output_{}.reshape({}, {}),'
        out_vec = 'output_{},'

        for i in range(self.num_arguments):
            lines['header_args'].append(hd.format('input', i))
            lines['python_args'].append(py.format('input', i))
            lines['c_args'].append(c.format('input', i))

        for i, matrix in enumerate(self.matrices):
            lines['header_args'].append(hd.format('output', i))
            lines['python_args'].append(py.format('output', i))
            lines['c_args'].append(c.format('output', i))
            nr, nc = matrix.shape
            if nc == 1:
                lines['output'].append(out_vec.format(i))
            else:
                lines['output'].append(out.format(i, nr, nc))

        indents = {'header_args': 18,
                   'python_args': 9,
                   'c_args': 13,
                   'output': 12}

        self.code_blocks = {k: wrap_and_indent(v, indents[k])[:-1] for k, v
                            in lines.items()}

    def doprint(self):
        """Returns the text of the four source files needed to compile
        Cython wrapper that evaluates the provided SymPy matrices.

        Returns
        =======
        setup_py : string
            The text of the setup.py file used to compile the Cython
            extension.
        cython_source : string
            The text of the Cython pyx file which includes the wrapper
            function ``eval``.
        c_header : string
            The text of the C header file that exposes the evaluate
            function.
        c_source : string
            The text of the C source file containing the function that
            evaluates the matrices.

        """

        c_header, c_source = self.c_matrix_generator.doprint(
            prefix=self.prefix + '_c')

        filling = {'prefix': self.prefix}
        filling.update(self.code_blocks)

        cython_source = self._pyx_template.format(**filling)
        setup_py = self._setup_py_template.format(**filling)

        return setup_py, cython_source, c_header, c_source

    def write(self, path=None):
        """Writes the four source files needed to compile the Cython
        function to the current working directory.

        Parameters
        ==========
        path : string
            The absolute or relative path to an existing directory to place
            the files instead of the cwd.

        """

        if path is None:
            path = os.getcwd()

        self.c_matrix_generator.write(self.prefix + '_c')

        setup_py, pyx, c_header, c_source = self.doprint()

        with open(os.path.join(path, self.prefix + '_setup.py'), 'w') as f:
            f.write(setup_py)

        with open(os.path.join(path, self.prefix + '.pyx'), 'w') as f:
            f.write(pyx)

    def compile(self, tmp_dir=None, verbose=False):
        """Returns a function which evaluates the matrices.

        Parameters
        ==========
        tmp_dir : string
            The path to an existing or non-existing directory where all of
            the generated files will be stored.
        verbose : boolean
            If true the output of the completed compilation steps will be
            printed.

        """

        base_prefix = self.prefix

        if tmp_dir is None:
            codedir = tempfile.mkdtemp(".pydy_compile")
        else:
            codedir = os.path.abspath(tmp_dir)

        if not os.path.exists(codedir):
            os.makedirs(codedir)

        self.prefix = '{}_{}'.format(base_prefix,
                                     CythonMatrixGenerator._module_counter)

        workingdir = os.getcwd()
        os.chdir(codedir)

        try:
            sys.path.append(codedir)
            self.write()
            cmd = [sys.executable, self.prefix + '_setup.py', 'build_ext',
                   '--inplace']
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            if verbose:
                print(output.decode())
            cython_module = importlib.import_module(self.prefix)
        except:
            raise Exception('Failed to compile and import Cython module.')
        finally:
            sys.path.remove(codedir)
            CythonMatrixGenerator._module_counter += 1
            os.chdir(workingdir)
            if tmp_dir is None:
                # rmtree fails on Windows with permissions errors, so skip the
                # removal on Windows.
                try:
                    shutil.rmtree(codedir)
                except OSError:
                    pass

        self.prefix = base_prefix

        return getattr(cython_module, 'eval')
