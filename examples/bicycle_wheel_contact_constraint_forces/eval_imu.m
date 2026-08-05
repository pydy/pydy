function [output_1] = eval_imu(input_1, input_2, input_3, input_4)
% function [output_1] = eval_imu(input_1, input_2, input_3, input_4)
%
% input_1 : [q4(t), q5(t), q7(t)]
% input_2 : [u3(t), u4(t), u5(t), u6(t), u7(t)]
% input_3 : [u3p(t), u4p(t), u5p(t), u6p(t), u7p(t)]
% input_4 : [bx, by, bz, d1, d2, d3, ex, ey, ez, g, lambda, rr]

    pydy_0 = input_2(2);
    pydy_1 = sin(input_4(11));
    pydy_2 = input_1(2);
    pydy_3 = sin(pydy_2);
    pydy_4 = cos(input_4(11));
    pydy_5 = cos(pydy_2);
    pydy_6 = pydy_1.*pydy_3 + pydy_4.*pydy_5;
    pydy_7 = pydy_0.*pydy_6;
    pydy_8 = pydy_1.*pydy_5;
    pydy_9 = pydy_3.*pydy_4;
    pydy_10 = pydy_8 - pydy_9;
    pydy_11 = input_2(1);
    pydy_12 = input_1(1);
    pydy_13 = cos(pydy_12);
    pydy_14 = pydy_11.*pydy_13;
    pydy_15 = pydy_10.*pydy_14 + pydy_7;
    pydy_16 = input_2(3);
    pydy_17 = sin(pydy_12);
    pydy_18 = pydy_11.*pydy_17;
    pydy_19 = pydy_16 + pydy_18;
    pydy_20 = -pydy_8 + pydy_9;
    pydy_21 = pydy_14.*pydy_6;
    pydy_22 = pydy_0.*pydy_20 + pydy_21;
    pydy_23 = input_3(3);
    pydy_24 = input_3(1);
    pydy_25 = pydy_0.*pydy_14;
    pydy_26 = pydy_17.*pydy_24 + pydy_23 + pydy_25;
    pydy_27 = input_4(10).*pydy_13;
    pydy_28 = input_4(12).*(pydy_19 + input_2(4));
    pydy_29 = pydy_0.^2.*input_4(12) + pydy_18.*pydy_28;
    pydy_30 = -pydy_25.*input_4(12) - input_4(12).*(pydy_26 + input_3(4));
    pydy_31 = -input_4(1).*pydy_19 + input_4(2).*pydy_15;
    pydy_32 = input_3(2);
    pydy_33 = pydy_13.*pydy_24;
    pydy_34 = -pydy_14.*pydy_16.*pydy_20 + pydy_16.*pydy_7 - ...
    pydy_18.*pydy_7 + pydy_20.*pydy_32 + pydy_33.*pydy_6;
    pydy_35 = input_4(1).*pydy_22 - input_4(3).*pydy_15;
    pydy_36 = -input_4(2).*pydy_22 + input_4(3).*pydy_19;
    pydy_37 = pydy_0.*pydy_10;
    pydy_38 = pydy_10.*pydy_33 - pydy_16.*pydy_21 + pydy_16.*pydy_37 - ...
    pydy_18.*pydy_37 + pydy_32.*pydy_6;
    pydy_39 = -pydy_14.*pydy_28 + pydy_32.*input_4(12);
    pydy_40 = input_1(3);
    pydy_41 = sin(pydy_40);
    pydy_42 = pydy_16.*pydy_41;
    pydy_43 = cos(pydy_40);
    pydy_44 = pydy_0.*pydy_5;
    pydy_45 = pydy_43.*pydy_44;
    pydy_46 = pydy_17.*pydy_41;
    pydy_47 = pydy_13.*pydy_43;
    pydy_48 = pydy_3.*pydy_47;
    pydy_49 = pydy_46 - pydy_48;
    pydy_50 = pydy_11.*pydy_49 + pydy_42 + pydy_45;
    pydy_51 = pydy_16.*pydy_43;
    pydy_52 = pydy_41.*pydy_44;
    pydy_53 = pydy_17.*pydy_43;
    pydy_54 = pydy_13.*pydy_41;
    pydy_55 = pydy_3.*pydy_54;
    pydy_56 = pydy_53 + pydy_55;
    pydy_57 = pydy_11.*pydy_56 + pydy_51 - pydy_52;
    pydy_58 = input_2(5);
    pydy_59 = pydy_0.*pydy_3;
    pydy_60 = pydy_14.*pydy_5;
    pydy_61 = pydy_59 + pydy_60;
    pydy_62 = pydy_58 + pydy_61;
    pydy_63 = pydy_29.*pydy_3;
    pydy_64 = -input_4(4).*pydy_19.^2 - input_4(4).*pydy_61.^2;
    pydy_65 = pydy_30.*pydy_5;
    pydy_66 = pydy_57.^2;
    pydy_67 = pydy_14.*pydy_3;
    pydy_68 = pydy_16.*pydy_44 - pydy_16.*pydy_67 - pydy_18.*pydy_44 + ...
    pydy_3.*pydy_32 + pydy_33.*pydy_5;
    pydy_69 = pydy_68 + input_3(5);
    pydy_70 = -input_4(5).*pydy_50 + input_4(6).*pydy_62;
    pydy_71 = input_4(7).*pydy_62 - input_4(9).*pydy_50;
    pydy_72 = pydy_44 - pydy_67;
    pydy_73 = input_4(4).*pydy_72;
    pydy_74 = input_4(4).*pydy_68 + pydy_19.*pydy_73;
    pydy_75 = -input_4(7).*pydy_57 + input_4(8).*pydy_50;
    pydy_76 = pydy_0.*pydy_11;
    pydy_77 = pydy_43.*pydy_58;
    pydy_78 = pydy_41.*pydy_58;
    pydy_79 = pydy_32.*pydy_5;
    pydy_80 = pydy_23.*pydy_43 + pydy_24.*pydy_56 - pydy_41.*pydy_79 + ...
    pydy_42.*pydy_59;
    pydy_81 = -pydy_19.*pydy_78 + pydy_42.*pydy_60 - pydy_72.*pydy_77 + ...
    pydy_76.*(-pydy_3.*pydy_46 + pydy_47) + pydy_80;
    pydy_82 = pydy_13.*pydy_5;
    pydy_83 = pydy_11.*(pydy_0.*pydy_47 + pydy_42.*pydy_82 - ...
    pydy_46.*pydy_58 - pydy_46.*pydy_59 + pydy_48.*pydy_58) - ...
    pydy_42.*pydy_58 - pydy_45.*pydy_58 + pydy_80;
    pydy_84 = -input_4(8).*pydy_62 + input_4(9).*pydy_57;
    pydy_85 = pydy_23.*pydy_41 + pydy_24.*pydy_49 + pydy_43.*pydy_79 - ...
    pydy_51.*pydy_59;
    pydy_86 = pydy_19.*pydy_77 - pydy_51.*pydy_60 - pydy_72.*pydy_78 + ...
    pydy_76.*(pydy_3.*pydy_53 + pydy_54) + pydy_85;

    output_1 = [pydy_15; pydy_19; pydy_22; -input_4(2).*pydy_34 + ...
    input_4(3).*pydy_26 - pydy_10.*pydy_27 + pydy_10.*pydy_29 + ...
    pydy_19.*pydy_31 - pydy_22.*pydy_35 + pydy_30.*pydy_6; ...
    input_4(1).*pydy_34 - input_4(3).*pydy_38 - input_4(10).*pydy_17 - ...
    pydy_15.*pydy_31 + pydy_22.*pydy_36 + pydy_39; -input_4(1).*pydy_26 + ...
    input_4(2).*pydy_38 + pydy_15.*pydy_35 - pydy_19.*pydy_36 + ...
    pydy_20.*pydy_30 - pydy_27.*pydy_6 + pydy_29.*pydy_6; pydy_50; pydy_57; ...
    pydy_62; input_4(5).*pydy_83 - input_4(6).*pydy_66 - ...
    input_4(8).*pydy_69 + input_4(9).*pydy_81 - input_4(10).*pydy_49 + ...
    pydy_39.*pydy_41 + pydy_41.*pydy_74 - pydy_43.*pydy_63 + ...
    pydy_43.*pydy_64 + pydy_43.*pydy_65 + pydy_57.*pydy_75 - ...
    pydy_62.*pydy_70 - pydy_62.*pydy_71; input_4(5).*pydy_57.*pydy_62 - ...
    input_4(5).*(pydy_11.*(pydy_0.*pydy_54 - pydy_51.*pydy_82 + ...
    pydy_53.*pydy_58 + pydy_53.*pydy_59 + pydy_55.*pydy_58) + ...
    pydy_51.*pydy_58 - pydy_52.*pydy_58 + pydy_85) + ...
    input_4(6).*pydy_50.*pydy_57 + input_4(6).*pydy_69 + ...
    input_4(7).*pydy_69 - input_4(9).*pydy_86 - input_4(10).*pydy_56 + ...
    pydy_39.*pydy_43 + pydy_41.*pydy_63 - pydy_41.*pydy_64 - ...
    pydy_41.*pydy_65 + pydy_43.*pydy_74 - pydy_50.*pydy_75 + ...
    pydy_62.*pydy_84; -input_4(4).*pydy_26 - input_4(5).*pydy_66 - ...
    input_4(6).*pydy_83 - input_4(7).*pydy_81 + input_4(8).*pydy_86 - ...
    pydy_27.*pydy_5 + pydy_29.*pydy_5 + pydy_3.*pydy_30 + pydy_50.*pydy_70 ...
    + pydy_50.*pydy_71 - pydy_57.*pydy_84 + pydy_61.*pydy_73];

end
