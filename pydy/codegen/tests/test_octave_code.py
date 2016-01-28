from sympy import ordered

from ...utils import sympy_equal_to_or_newer_than
from ...models import multi_mass_spring_damper


if sympy_equal_to_or_newer_than('0.7.6'):

    from ..octave_code import OctaveMatrixGenerator

    def test_OctaveMatrixGenerator():

        expected_m_file = """\
function [output_1] = eval_mats(input_1, input_2, input_3)
% function [output_1] = eval_mats(input_1, input_2, input_3)
%
% input_1 : [x0(t), x1(t), x2(t)]
% input_2 : [v0(t), v1(t), v2(t)]
% input_3 : [c0, c1, c2, k0, k1, k2, m0, m1, m2]

    pydy_0 = input_2(1);
    pydy_1 = input_2(2);
    pydy_2 = input_2(3);
    pydy_3 = input_3(8) + input_3(9);
    pydy_4 = 1./(input_3(7) + pydy_3);
    pydy_5 = -input_3(1).*pydy_0 - input_3(4).*input_1(1);
    pydy_6 = input_3(9).*pydy_4;
    pydy_7 = input_3(9) - pydy_3.*pydy_6;
    pydy_8 = 1./(-pydy_3.^2.*pydy_4 + pydy_3);
    pydy_9 = -input_3(2).*pydy_1 - input_3(5).*input_1(2) - ...
    pydy_3.*pydy_4.*pydy_5;
    pydy_10 = (-input_3(3).*pydy_2 - input_3(6).*input_1(3) - ...
    pydy_5.*pydy_6 - pydy_7.*pydy_8.*pydy_9)./(-input_3(9).^2.*pydy_4 + ...
    input_3(9) - pydy_7.^2.*pydy_8);
    pydy_11 = pydy_8.*(-pydy_10.*pydy_7 + pydy_9);

    output_1 = [pydy_0; pydy_1; pydy_2; pydy_4.*(-input_3(9).*pydy_10 - ...
    pydy_11.*pydy_3 + pydy_5); pydy_11; pydy_10];

end
"""

        sys = multi_mass_spring_damper(3)
        q = sys.coordinates
        u = sys.speeds
        p = list(ordered(sys.constants_symbols))
        sym_rhs = sys.eom_method.rhs()
        g = OctaveMatrixGenerator([q, u, p], [sym_rhs])
        assert g.doprint() == expected_m_file
