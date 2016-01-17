from pydy.models import multi_mass_spring_damper
from pydy.codegen.matrix import OctaveMatrixGenerator


def test_OctaveMatrixGenerator():

    expected_m_file = """\
function [output_1] = eval_mats(input_1, ...
                                input_2, ...
                                input_3)

    pydy_0 = input_2(1);
    pydy_1 = input_2(2);
    pydy_2 = input_2(3);
    pydy_3 = input_3(3) + input_3(8);
    pydy_4 = 1./(input_3(6) + pydy_3);
    pydy_5 = -input_3(4).*pydy_0 - input_3(9).*input_1(1);
    pydy_6 = input_3(8).*pydy_4;
    pydy_7 = input_3(8) - pydy_3.*pydy_6;
    pydy_8 = 1./(-pydy_3.^2.*pydy_4 + pydy_3);
    pydy_9 = -input_3(1).*pydy_1 - input_3(5).*input_1(2) - ...
    pydy_3.*pydy_4.*pydy_5;
    pydy_10 = (-input_3(7).*pydy_2 - input_3(2).*input_1(3) - pydy_5.*pydy_6 - ...
    pydy_7.*pydy_8.*pydy_9)./(-input_3(8).^2.*pydy_4 + input_3(8) - ...
    pydy_7.^2.*pydy_8);
    pydy_11 = pydy_8.*(-pydy_10.*pydy_7 + pydy_9);

    output_1 = [pydy_0; pydy_1; pydy_2; pydy_4.*(-input_3(8).*pydy_10 -
    pydy_11.*pydy_3 + pydy_5); pydy_11; pydy_10];
end"""

    sys = multi_mass_spring_damper(3)
    q = sys.coordinates
    u = sys.speeds
    p = sys.constants_symbols
    sym_rhs = sys.eom_method.rhs()
    g = OctaveMatrixGenerator([q, u, p], [sym_rhs])
    assert g.doprint() == expected_m_file
