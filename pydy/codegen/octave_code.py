#!/usr/bin/env python

import os

try:
    from sympy.printing.octave import OctaveCodePrinter
except ImportError:
    raise ImportError('Octave code printing is only available with SymPy >= 0.7.6.')

from .matrix_generator import MatrixGenerator


class OctaveMatrixGenerator(MatrixGenerator):
    """This class generates Octave/Matlab source files that simultaneously
    numerically evaluate any number of SymPy matrices.

    """

    _idx_start = 1
    _idx_delim = "()"
    _base_printer = OctaveCodePrinter
    _type_declar = ''  # prepended to variable introductions
    _line_contin = ' ...'
    _comment_char = '%'

    # TODO : The first two lines will not wrap. For many inputs/outputs it
    # would be nice to have some wrapping.
    _m_template = """\
function [{output_args}] = {prefix}({input_args})
% function [{output_args}] = {prefix}({input_args})
%
{docstring}

{subexprs}

{outputs}

end
"""

    def doprint(self, prefix='eval_mats'):
        """Returns a string that implements the function.

        Parameters
        ==========
        prefix : string, optional
            The name of the Octave/Matlab function.

        """
        self.code_blocks['prefix'] = prefix

        return self._m_template.format(**self.code_blocks)

    def write(self, prefix='eval_mats', path=None):
        """Writes the <prefix>.m file to disc at the give path location."""

        if path is None:
            path = os.getcwd()

        text = self.doprint(prefix=prefix)

        with open(os.path.join(path, prefix + '.m'), 'w') as f:
            f.write(text)
