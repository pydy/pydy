#!/usr/bin/env python

"""This module contains source code dedicated to generating C code from
matrices generated from sympy.physics.mechanics."""

import os
import itertools

import sympy as sm
import sympy.physics.mechanics as me
from sympy.printing.ccode import CCodePrinter

from ..utils import wrap_and_indent, find_dynamicsymbols


class CMatrixGenerator(object):
    """This class generates C source files that simultaneously numerically
    evaluate any number of SymPy matrices.

    """

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

    def __init__(self, arguments, matrices):
        """

        Parameters
        ==========
        arguments : sequences of sequences of SymPy Symbol or Function.
            Each of the sequences will be converted to input arrays in the C
            function. All of the symbols/functions contained in ``matrices``
            need to be in the sequences, but the sequences can also contain
            extra symbols/functions that are not contained in the matrices.
        matrices : sequence of SymPy.Matrix
            A sequence of the matrices that should be evaluated in the
            function. The expressions should contain only sympy.Symbol or
            sympy.Function that are functions of me.dynamicsymbols._t.

        """

        required_args = set()

        for matrix in matrices:
            # TODO : SymPy 0.7.4 does not have Matrix.free_symbols so we
            # manually compute them instead of calling:
            # required_args |= matrix.free_symbols
            required_args |= set().union(*[i.free_symbols for i in matrix])
            required_args |= find_dynamicsymbols(matrix)

        required_args.remove(me.dynamicsymbols._t)

        all_arguments = set(itertools.chain(*arguments))

        for required_arg in required_args:
            if required_arg not in all_arguments:
                msg = "{} is missing from the argument sequences."
                raise ValueError(msg.format(required_arg))

        self.matrices = matrices
        self.arguments = arguments

        self._generate_cse()
        self._generate_code_blocks()

    def _generate_cse(self):

        # This makes a big long list of every expression in all of the
        # matrices.
        exprs = []
        for matrix in self.matrices:
            exprs += matrix[:]

        # Compute the common subexpresions.
        gen = sm.numbered_symbols('pydy_')
        self.subexprs, simplified_exprs = sm.cse(exprs, symbols=gen)

        # Turn the expressions back into matrices of the same type.
        simplified_matrices = []
        idx = 0
        for matrix in self.matrices:
            num_rows, num_cols = matrix.shape
            length = num_rows * num_cols
            m = type(matrix)(simplified_exprs[idx:idx + length])
            simplified_matrices.append(m.reshape(num_rows, num_cols))
            idx += length

        self.simplified_matrices = tuple(simplified_matrices)

    def _generate_pydy_c_printer(self):
        """Returns a subclass of sympy.printing.CCodePrinter to print
        appropriate C array index calls for all of the symbols in the
        equations of motion.

        Examples
        --------

        >>> from sympy import symbols
        >>> from sympy.physics.mechanics import dynamicsymbols
        >>> from pydy.codegen.c_code import CMatrixGenerator
        >>> generator = CMatrixGenerator(...)
        >>> PyDyCCodePrinter = generator._generate_pydy_c_printer()
        >>> printer = PyDyCCodePrinter()
        >>> m = symbols('m') # m is the first constant in the EoMs
        >>> printer.doprint(m)
        input_0[0]
        >>> q = dynamicsymbols('q') # q is the second coordinate in the EoMs
        >>> printer.doprint(q)
        input_1[1]
        >>> F = dynamicsymbols('F') # F is the third specified in the EoMs
        >>> printer.doprint(F)
        input_3[2]

        """

        array_index_map = {}
        for i, arg_set in enumerate(self.arguments):
            for j, var in enumerate(arg_set):
                array_index_map[var] = r'input_{}[{}]'.format(i, j)

        class PyDyCCodePrinter(CCodePrinter):

            def _print_Function(self, e):
                if e in array_index_map.keys():
                    return array_index_map[e]
                else:
                    return super(PyDyCCodePrinter, self)._print_Function(e)

            def _print_Symbol(self, e):
                if e in array_index_map.keys():
                    return array_index_map[e]
                else:
                    return super(PyDyCCodePrinter, self)._print_Symbol(e)

        return PyDyCCodePrinter

    def comma_lists(self):
        """Returns a string output for each of the sequences of SymPy
        arguments."""

        comma_lists = []

        for i, arg_set in enumerate(self.arguments):
            comma_lists.append(', '.join([str(s) for s in arg_set]))

        return tuple(comma_lists)

    def _generate_code_blocks(self):
        """Writes the blocks of code for the C file."""

        printer = self._generate_pydy_c_printer()()

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
