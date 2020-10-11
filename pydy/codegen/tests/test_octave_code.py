from pkg_resources import parse_version
from sympy import ordered, __version__

from ...models import multi_mass_spring_damper
from ..octave_code import OctaveMatrixGenerator

SYMPY_VERSION = __version__
del __version__


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

# SymPy > 1.0 outputs different cse results.
    expected_m_file_new = """\
function [output_1] = eval_mats(input_1, input_2, input_3)
% function [output_1] = eval_mats(input_1, input_2, input_3)
%
% input_1 : [x0(t), x1(t), x2(t)]
% input_2 : [v0(t), v1(t), v2(t)]
% input_3 : [c0, c1, c2, k0, k1, k2, m0, m1, m2]

    pydy_0 = input_2(1);
    pydy_1 = input_2(2);
    pydy_2 = input_2(3);
    pydy_3 = 1./(input_3(7) + input_3(8) + input_3(9));
    pydy_4 = -input_3(1).*pydy_0;
    pydy_5 = -input_3(4).*input_1(1);
    pydy_6 = input_3(8) + input_3(9);
    pydy_7 = -input_3(9).*pydy_3.*pydy_6 + input_3(9);
    pydy_8 = 1./(input_3(8) + input_3(9) - pydy_3.*pydy_6.^2);
    pydy_9 = 1./(-input_3(9).^2.*pydy_3 + input_3(9) - pydy_7.^2.*pydy_8);
    pydy_10 = pydy_4 + pydy_5;
    pydy_11 = -input_3(2).*pydy_1;
    pydy_12 = -input_3(5).*input_1(2);
    pydy_13 = -pydy_10.*pydy_3.*pydy_6;
    pydy_14 = -input_3(3).*pydy_2 - input_3(6).*input_1(3) - ...
    input_3(9).*pydy_10.*pydy_3 - pydy_7.*pydy_8.*(pydy_11 + pydy_12 + ...
    pydy_13);
    pydy_15 = pydy_11 + pydy_12 + pydy_13 - pydy_14.*pydy_7.*pydy_9;

    output_1 = [pydy_0; pydy_1; pydy_2; ...
    pydy_3.*(-input_3(9).*pydy_14.*pydy_9 - pydy_15.*pydy_6.*pydy_8 + ...
    pydy_4 + pydy_5); pydy_15.*pydy_8; pydy_14.*pydy_9];

end
"""

    expected_m_file_1p1 = """\
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

    expected_m_file_1p2 = """\
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
    pydy_5 = input_3(9).*pydy_4;
    pydy_6 = input_3(9) - pydy_3.*pydy_5;
    pydy_7 = 1./(-pydy_3.^2.*pydy_4 + pydy_3);
    pydy_8 = -input_3(1).*pydy_0 - input_3(4).*input_1(1);
    pydy_9 = -input_3(2).*pydy_1 - input_3(5).*input_1(2) - ...
    pydy_3.*pydy_4.*pydy_8;
    pydy_10 = (-input_3(3).*pydy_2 - input_3(6).*input_1(3) - ...
    pydy_5.*pydy_8 - pydy_6.*pydy_7.*pydy_9)./(-input_3(9).^2.*pydy_4 + ...
    input_3(9) - pydy_6.^2.*pydy_7);
    pydy_11 = pydy_7.*(-pydy_10.*pydy_6 + pydy_9);

    output_1 = [pydy_0; pydy_1; pydy_2; pydy_4.*(-input_3(9).*pydy_10 - ...
    pydy_11.*pydy_3 + pydy_8); pydy_11; pydy_10];

end
"""

    sys = multi_mass_spring_damper(3)
    q = sys.coordinates
    u = sys.speeds
    p = list(ordered(sys.constants_symbols))
    sym_rhs = sys.eom_method.rhs()
    g = OctaveMatrixGenerator([q, u, p], [sym_rhs])
    if parse_version(SYMPY_VERSION) >= parse_version('1.2'):
        assert g.doprint() == expected_m_file_1p2
    elif parse_version(SYMPY_VERSION) >= parse_version('1.1'):
        assert g.doprint() == expected_m_file_1p1
    elif parse_version(SYMPY_VERSION) > parse_version('1.0'):
        assert g.doprint() == expected_m_file_new
    else:
        assert g.doprint() == expected_m_file
