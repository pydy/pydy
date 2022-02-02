#!/usr/bin/env python

"""This module contains source code dedicated to generating C code from
matrices generated from sympy.physics.mechanics."""

import os

from pkg_resources import parse_version
import sympy as sm

SYMPY_VERSION = sm.__version__

try:
    try:
        if parse_version(SYMPY_VERSION) >= parse_version('1.7'):
            from sympy.printing.c import C99CodePrinter as CCodePrinter
        else:
            from sympy.printing.ccode import C99CodePrinter as CCodePrinter
    except ImportError:
        # SymPy 1.0 and lower uses this version.
        from sympy.printing.ccode import CCodePrinter
except ModuleNotFoundError:
    # SymPy > 1.5 renamed the module ccode to c
    from sympy.printing.c import C99CodePrinter as CCodePrinter

from .matrix_generator import MatrixGenerator
from ..utils import wrap_and_indent


class CMatrixGenerator(MatrixGenerator):
    """This class generates C source files that simultaneously numerically
    evaluate any number of SymPy matrices.

    """

    _base_printer = CCodePrinter

    _c_header_template = """\
void evaluate(
{input_args}
{output_args}
             );
/*

{input_docstring}

*/"""

    _c_source_template = """\
#include <math.h>{header_include}

void evaluate(
{input_args}
{output_args}
             )
{{

{subexprs}

{outputs}

}}\
"""

    def _generate_code_blocks(self):
        """Writes the blocks of code for the C file."""

        # TODO : This could use the super classes method with some tweaks.

        printer = self._generate_pydy_printer()()

        self.code_blocks = {}

        lines = []
        for i, input_arg in enumerate(self.arguments):
            lines.append('double input_{}[{}],'.format(i, len(input_arg)))
        self.code_blocks['input_args'] = wrap_and_indent(lines, 14)

        lines = []
        for i, output_arg in enumerate(self.matrices):
            nr, nc = output_arg.shape
            lines.append('double output_{}[{}],'.format(i, nr * nc))
            self.code_blocks['output_args'] = \
                wrap_and_indent(lines, 14)[:-1]  # remove last comma

        lines = []
        for i, (input_arg, explan) in enumerate(zip(self.arguments,
                                                    self.comma_lists())):
            lines.append('input_{}[{}] : [{}]'.format(i, len(input_arg),
                                                      explan))
        self.code_blocks['input_docstring'] = wrap_and_indent(lines, 0)

        lines = []
        for var, expr in self.subexprs:
            var_str = printer.doprint(var)
            expr_str = printer.doprint(expr)
            lines.append('double {} = {};'.format(var_str, expr_str))
        self.code_blocks['subexprs'] = wrap_and_indent(lines)

        outputs = ''
        for i, output in enumerate(self.simplified_matrices):
            nr, nc = output.shape
            lhs = sm.MatrixSymbol('output_{}'.format(i), nr, nc)
            try:
                code_str = printer.doprint(output, lhs)
            except AttributeError:
                # The above fails in SymPy 0.7.4.1 because Matrix printing
                # isn't supported.
                code_lines = []
                for j, element in enumerate(output):
                    assignment = 'output_{}[{}]'.format(i, j)
                    code_lines.append(printer.doprint(element, assignment))
                code_str = '\n'.join(code_lines)
            outputs += wrap_and_indent(code_str.split('\n'))
            if i != len(self.simplified_matrices) - 1:
                outputs += '\n\n'  # space between each output

        self.code_blocks['outputs'] = outputs

    def doprint(self, prefix=None):
        """Returns a string each for the header and the source files.

        Parameters
        ==========
        prefix : string, optional
            A prefix for the name of the header file. This will cause an
            include statement to to be added to the source.

        """

        if prefix is not None:
            filling = {'header_include': '\n#include "{}.h"'.format(prefix)}
        else:
            filling = {'header_include': ''}

        filling.update(self.code_blocks)

        c_header = self._c_header_template.format(**filling)
        c_source = self._c_source_template.format(**filling)

        return c_header, c_source

    def write(self, prefix, path=None):
        """Writes a header and source file to disk.

        Parameters
        ==========
        prefix : string
            Two files will be generated: ``<prefix>.c`` and
            ``<prefix>.h``.

        """

        if path is None:
            path = os.getcwd()

        header, source = self.doprint(prefix=prefix)

        with open(os.path.join(path, prefix + '.h'), 'w') as f:
            f.write(header)

        with open(os.path.join(path, prefix + '.c'), 'w') as f:
            f.write(source)


class _CLUsolveGenerator(CMatrixGenerator):
    """This is a private undocumented class that supports the
    ``linear_sys_solver='sympy'`` in CythonMatrixGenerator. It cse's A and b of
    a linear system Ax=b, then solves the linear system symbolically and cse's
    the result x. This is a more efficient way to get the symbolic solution of
    a linear system encoded in generated C code."""

    def _generate_cse(self, prefix='pydy_'):
        # NOTE : This assumes the first two items in self.matrices are A and b
        # of and Ax=b system. This also ignores cse=False.

        gen1 = sm.numbered_symbols(prefix)
        subexprs1, mats_simp = sm.cse(self.matrices, symbols=gen1)

        A_simp = mats_simp[0]
        b_simp = mats_simp[1]

        x = A_simp.LUsolve(b_simp)

        gen2 = sm.numbered_symbols(prefix, start=len(subexprs1))
        subexprs2, x_simp = sm.cse(x, symbols=gen2)

        # swap the b matrix with the x result
        mats_simp[1] = x_simp[0]

        self.subexprs = subexprs1 + subexprs2
        self.simplified_matrices = tuple(mats_simp)
