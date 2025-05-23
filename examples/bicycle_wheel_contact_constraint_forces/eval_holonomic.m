function [output_1] = eval_holonomic(input_1, input_2)
% function [output_1] = eval_holonomic(input_1, input_2)
%
% input_1 : [q4(t), q5(t), q7(t)]
% input_2 : [d1, d2, d3, rf, rr]

    pydy_0 = input_1(1);
    pydy_1 = cos(pydy_0);
    pydy_2 = input_1(2);
    pydy_3 = sin(pydy_2);
    pydy_4 = sin(pydy_0);
    pydy_5 = input_1(3);
    pydy_6 = sin(pydy_5);
    pydy_7 = pydy_4.*pydy_6;
    pydy_8 = cos(pydy_5);
    pydy_9 = pydy_1.*pydy_8;
    pydy_10 = cos(pydy_2);
    pydy_11 = -pydy_3.*pydy_7 + pydy_9;
    pydy_12 = pydy_10.^2.*pydy_7 - pydy_11.*pydy_3;
    pydy_13 = pydy_10.*pydy_11.*pydy_8 + pydy_10.*pydy_6.*(pydy_1.*pydy_6 + ...
    pydy_3.*pydy_4.*pydy_8);
    pydy_14 = input_2(4)./sqrt(pydy_12.^2 + pydy_13.^2);

    output_1 = -input_2(1).*pydy_1.*pydy_3 + pydy_1.*pydy_10.*(input_2(2) + ...
    pydy_13.*pydy_14) - pydy_1.*input_2(5) + (input_2(3) + ...
    pydy_12.*pydy_14).*(-pydy_3.*pydy_9 + pydy_7);

end
