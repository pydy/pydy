function [output_1, output_2] = eval_dep_speeds(input_1, input_2, input_3)
% function [output_1, output_2] = eval_dep_speeds(input_1, input_2, input_3)
%
% input_1 : [q4(t), q5(t), q7(t)]
% input_2 : [u4(t), u6(t), u7(t)]
% input_3 : [d1, d2, d3, rf, rr]

    pydy_0 = input_1(1);
    pydy_1 = sin(pydy_0);
    pydy_2 = pydy_1.*input_3(5);
    pydy_3 = input_1(2);
    pydy_4 = sin(pydy_3);
    pydy_5 = input_3(1).*pydy_4;
    pydy_6 = pydy_1.*pydy_5;
    pydy_7 = input_1(3);
    pydy_8 = cos(pydy_7);
    pydy_9 = pydy_1.*pydy_8;
    pydy_10 = sin(pydy_7);
    pydy_11 = cos(pydy_0);
    pydy_12 = pydy_10.*pydy_11;
    pydy_13 = pydy_12.*pydy_4 + pydy_9;
    pydy_14 = input_3(3).*pydy_13;
    pydy_15 = cos(pydy_3);
    pydy_16 = pydy_1.*pydy_10;
    pydy_17 = pydy_11.*pydy_8;
    pydy_18 = -pydy_16.*pydy_4 + pydy_17;
    pydy_19 = pydy_15.^2.*pydy_16 - pydy_18.*pydy_4;
    pydy_20 = pydy_12 + pydy_4.*pydy_9;
    pydy_21 = pydy_10.*pydy_15;
    pydy_22 = pydy_15.*pydy_8;
    pydy_23 = pydy_18.*pydy_22 + pydy_20.*pydy_21;
    pydy_24 = input_3(4)./sqrt(pydy_19.^2 + pydy_23.^2);
    pydy_25 = pydy_13.*pydy_24;
    pydy_26 = pydy_19.*pydy_25;
    pydy_27 = -pydy_14 - pydy_26;
    pydy_28 = input_3(2).*pydy_13 + pydy_23.*pydy_25;
    pydy_29 = pydy_16 - pydy_17.*pydy_4;
    pydy_30 = input_3(2).*pydy_29;
    pydy_31 = pydy_11.*pydy_15;
    pydy_32 = input_3(3).*pydy_31;
    pydy_33 = pydy_19.*pydy_24;
    pydy_34 = pydy_31.*pydy_33;
    pydy_35 = pydy_23.*pydy_24;
    pydy_36 = pydy_29.*pydy_35;
    pydy_37 = -pydy_34 + pydy_36;
    pydy_38 = input_3(3).*pydy_8;
    pydy_39 = pydy_33.*pydy_8;
    pydy_40 = -pydy_38 - pydy_39;
    pydy_41 = input_3(2).*pydy_10;
    pydy_42 = pydy_10.*pydy_35;
    pydy_43 = input_3(2).*pydy_8;
    pydy_44 = pydy_35.*pydy_8 + pydy_43;
    pydy_45 = pydy_22.*pydy_35;
    pydy_46 = input_3(1).*pydy_15;
    pydy_47 = pydy_1.*pydy_15;
    pydy_48 = -pydy_30 + pydy_32 + pydy_34 - pydy_36;
    pydy_49 = -pydy_41 - pydy_42;
    pydy_50 = input_2(1);
    pydy_51 = pydy_21.*pydy_50;
    pydy_52 = input_3(3).*pydy_51 + pydy_33.*pydy_51;
    pydy_53 = pydy_15.*pydy_50;
    pydy_54 = -pydy_35.*pydy_51 - pydy_41.*pydy_53;
    pydy_55 = pydy_4.*pydy_50 + input_2(3);
    pydy_56 = input_3(3).*pydy_55 + pydy_33.*pydy_55 - pydy_43.*pydy_53 - ...
    pydy_45.*pydy_50;
    pydy_57 = pydy_11.*pydy_50;

    output_1 = [-pydy_2 + pydy_21.*(pydy_30 - pydy_32 + pydy_37) + ...
    pydy_22.*pydy_28 + pydy_27.*pydy_4 - pydy_6 pydy_21.*(pydy_41 + ...
    pydy_42) + pydy_22.*pydy_44 + pydy_4.*pydy_40 - pydy_5 - input_3(5) ...
    -pydy_33.*pydy_4 + pydy_45; pydy_1.^2.*pydy_46 + pydy_11.^2.*pydy_46 + ...
    pydy_18.*pydy_48 + pydy_20.*pydy_28 + pydy_47.*(pydy_14 + pydy_26) ...
    pydy_1.*pydy_46 + pydy_18.*pydy_49 + pydy_20.*pydy_44 + ...
    pydy_47.*(pydy_38 + pydy_39) pydy_20.*pydy_35 + pydy_33.*pydy_47; ...
    pydy_13.*pydy_48 + pydy_27.*pydy_31 + pydy_28.*pydy_29 ...
    -input_3(1).*pydy_31 + pydy_13.*pydy_49 + pydy_29.*pydy_44 + ...
    pydy_31.*pydy_40 pydy_37];

    output_2 = [pydy_21.*pydy_56 - pydy_22.*pydy_54 - pydy_4.*pydy_52 + ...
    input_3(5).*input_2(2); -pydy_18.*pydy_56 - pydy_20.*pydy_54 + ...
    pydy_47.*pydy_52 - pydy_5.*pydy_57 - pydy_57.*input_3(5); ...
    -pydy_13.*pydy_56 - pydy_2.*pydy_50 - pydy_29.*pydy_54 - ...
    pydy_31.*pydy_52 - pydy_50.*pydy_6];

end
