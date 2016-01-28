#!/usr/bin/env python

import itertools

import sympy as sm
import sympy.physics.mechanics as me
from sympy.printing.codeprinter import CodePrinter

from ..utils import wrap_and_indent, find_dynamicsymbols


class MatrixGenerator(object):
    """This abstract base class generates source files that simultaneously
    numerically evaluate any number of SymPy matrices.

    """

    _idx_start = 0
    _idx_delim = "[]"
    _base_printer = CodePrinter
    # TODO : I opened an issue with SymPy to deal with the type declaration,
    # https://github.com/sympy/sympy/issues/10446.
    _type_declar = 'double '  # prepended to variable introductions
    _line_contin = None
    _comment_char = '#'

    def __init__(self, arguments, matrices, cse=True):
        """

        Parameters
        ==========
        arguments : sequence of sequences of SymPy Symbol or Function
            Each of the sequences will be converted to input arrays in the
            generated function. All of the symbols/functions contained in
            ``matrices`` need to be in the sequences, but the sequences can
            also contain extra symbols/functions that are not contained in the
            matrices.
        matrices : sequence of SymPy.Matrix
            A sequence of the matrices that should be evaluated in the
            function. The expressions should contain only sympy.Symbol or
            sympy.Function that are functions of me.dynamicsymbols._t.
        cse : boolean
            Find and replace common sub-expressions in ``matrices`` if True.

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

        if cse:
            self._generate_cse()
        else:
            self._ignore_cse()

        self._generate_code_blocks()

    def _generate_cse(self, prefix='pydy_'):

        # This makes a big long list of every expression in all of the
        # matrices.
        exprs = []
        for matrix in self.matrices:
            exprs += matrix[:]

        # Compute the common subexpresions.
        gen = sm.numbered_symbols(prefix)
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

    def _ignore_cse(self):
        self.subexprs = []
        self.simplified_matrices = self.matrices

    def _generate_pydy_printer(self):
        """Returns a subclass of a sympy.printing.Printer to print appropriate
        array index calls for all of the symbols in the equations of motion."""

        array_index_map = {}
        for i, arg_set in enumerate(self.arguments):
            for j, var in enumerate(arg_set):
                array_index_map[var] = r'input_{}{}{}{}'.format(
                    i + self._idx_start, self._idx_delim[0],
                    j + self._idx_start, self._idx_delim[1])

        class PyDyCodePrinter(self._base_printer):

            def _print_Function(self, e):
                if e in array_index_map.keys():
                    return array_index_map[e]
                else:
                    return super(PyDyCodePrinter, self)._print_Function(e)

            def _print_Symbol(self, e):
                if e in array_index_map.keys():
                    return array_index_map[e]
                else:
                    return super(PyDyCodePrinter, self)._print_Symbol(e)

        return PyDyCodePrinter

    def comma_lists(self):
        """Returns a string output for each of the sequences of SymPy
        arguments."""

        comma_lists = []

        for i, arg_set in enumerate(self.arguments):
            comma_lists.append(', '.join([str(s) for s in arg_set]))

        return tuple(comma_lists)

    def _generate_code_blocks(self):
        """Writes the blocks of code for the source file."""

        printer = self._generate_pydy_printer()()

        self.code_blocks = {}

        # Creates a string of comma separated inputs.
        lines = []
        for i, input_arg in enumerate(self.arguments):
            lines.append('{}input_{}'.format(self._type_declar,
                                             i + self._idx_start))
        self.code_blocks['input_args'] = ', '.join(lines)

        # Creates a string of comma separated outputs.
        lines = []
        for i, output_arg in enumerate(self.matrices):
            lines.append('{}output_{}'.format(self._type_declar,
                                              i + self._idx_start))
        self.code_blocks['output_args'] = ', '.join(lines)

        lines = []
        for i, (input_arg, explan) in enumerate(zip(self.arguments,
                                                    self.comma_lists())):
            lines.append('{} input_{} : [{}]'.format(self._comment_char, i +
                                                     self._idx_start, explan))
        self.code_blocks['docstring'] = wrap_and_indent(lines, 0,
                                                        comment=self._comment_char)

        lines = []
        for var, expr in self.subexprs:
            lines.append(printer.doprint(expr, var))
        self.code_blocks['subexprs'] = wrap_and_indent(lines,
                                                       continuation=self._line_contin)

        outputs = ''
        for i, output in enumerate(self.simplified_matrices):
            nr, nc = output.shape
            lhs = sm.MatrixSymbol('output_{}'.format(i + self._idx_start), nr, nc)
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
            outputs += wrap_and_indent(code_str.split('\n'),
                                       continuation=self._line_contin)
            if i != len(self.simplified_matrices) - 1:
                outputs += '\n\n'  # space between each output

        self.code_blocks['outputs'] = outputs
