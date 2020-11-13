function [output_1] = eval_imu(input_1, input_2, input_3, input_4)
% function [output_1] = eval_imu(input_1, input_2, input_3, input_4)
%
% input_1 : [q4(t), q5(t), q7(t)]
% input_2 : [u3(t), u4(t), u5(t), u6(t), u7(t)]
% input_3 : [u3p(t), u4p(t), u5p(t), u6p(t), u7p(t)]
% input_4 : [bx, bz, d1, d2, d3, ex, ez, g, lambda, rr]

    pydy_0 = input_2(2);
    pydy_1 = sin(input_4(9));
    pydy_2 = input_1(2);
    pydy_3 = sin(pydy_2);
    pydy_4 = cos(input_4(9));
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
    pydy_23 = pydy_19.^2;
    pydy_24 = input_3(3);
    pydy_25 = input_3(1);
    pydy_26 = pydy_0.*pydy_14;
    pydy_27 = pydy_17.*pydy_25 + pydy_24 + pydy_26;
    pydy_28 = pydy_10.*pydy_13;
    pydy_29 = input_4(10).*(pydy_19 + input_2(4));
    pydy_30 = pydy_0.^2.*input_4(10) + pydy_18.*pydy_29;
    pydy_31 = -pydy_26.*input_4(10) - input_4(10).*(pydy_27 + input_3(4));
    pydy_32 = input_4(1).*pydy_22 - input_4(2).*pydy_15;
    pydy_33 = input_3(2);
    pydy_34 = pydy_13.*pydy_6;
    pydy_35 = pydy_0.*pydy_10;
    pydy_36 = -pydy_14.*pydy_29 + pydy_33.*input_4(10);
    pydy_37 = input_1(3);
    pydy_38 = sin(pydy_37);
    pydy_39 = pydy_16.*pydy_38;
    pydy_40 = cos(pydy_37);
    pydy_41 = pydy_0.*pydy_5;
    pydy_42 = pydy_40.*pydy_41;
    pydy_43 = pydy_17.*pydy_38;
    pydy_44 = pydy_13.*pydy_40;
    pydy_45 = pydy_3.*pydy_44;
    pydy_46 = pydy_43 - pydy_45;
    pydy_47 = pydy_11.*pydy_46 + pydy_39 + pydy_42;
    pydy_48 = pydy_16.*pydy_40;
    pydy_49 = pydy_38.*pydy_41;
    pydy_50 = pydy_17.*pydy_40;
    pydy_51 = pydy_13.*pydy_38;
    pydy_52 = pydy_3.*pydy_51;
    pydy_53 = pydy_50 + pydy_52;
    pydy_54 = pydy_11.*pydy_53 + pydy_48 - pydy_49;
    pydy_55 = input_2(5);
    pydy_56 = pydy_0.*pydy_3;
    pydy_57 = pydy_14.*pydy_5;
    pydy_58 = pydy_56 + pydy_57;
    pydy_59 = pydy_55 + pydy_58;
    pydy_60 = pydy_3.*pydy_30;
    pydy_61 = -input_4(3).*pydy_23 - input_4(3).*pydy_58.^2;
    pydy_62 = pydy_31.*pydy_5;
    pydy_63 = pydy_54.^2;
    pydy_64 = -input_4(4).*pydy_47 + input_4(5).*pydy_59;
    pydy_65 = input_4(6).*pydy_59 - input_4(7).*pydy_47;
    pydy_66 = pydy_14.*pydy_3;
    pydy_67 = pydy_41 - pydy_66;
    pydy_68 = input_4(3).*pydy_67;
    pydy_69 = pydy_13.*pydy_5;
    pydy_70 = pydy_16.*pydy_41 - pydy_16.*pydy_66 - pydy_18.*pydy_41 + ...
    pydy_25.*pydy_69 + pydy_3.*pydy_33;
    pydy_71 = input_4(3).*pydy_70 + pydy_19.*pydy_68;
    pydy_72 = pydy_0.*pydy_11;
    pydy_73 = pydy_40.*pydy_55;
    pydy_74 = pydy_38.*pydy_55;
    pydy_75 = pydy_33.*pydy_5;
    pydy_76 = pydy_24.*pydy_40 + pydy_25.*pydy_53 - pydy_38.*pydy_75 + ...
    pydy_39.*pydy_56;
    pydy_77 = -pydy_19.*pydy_74 + pydy_39.*pydy_57 - pydy_67.*pydy_73 + ...
    pydy_72.*(-pydy_3.*pydy_43 + pydy_44) + pydy_76;
    pydy_78 = pydy_11.*(pydy_0.*pydy_44 + pydy_39.*pydy_69 - ...
    pydy_43.*pydy_55 - pydy_43.*pydy_56 + pydy_45.*pydy_55) - ...
    pydy_39.*pydy_55 - pydy_42.*pydy_55 + pydy_76;
    pydy_79 = pydy_70 + input_3(5);
    pydy_80 = pydy_54.*pydy_59;
    pydy_81 = pydy_47.*pydy_54;
    pydy_82 = pydy_24.*pydy_38 + pydy_25.*pydy_46 + pydy_40.*pydy_75 - ...
    pydy_48.*pydy_56;

    output_1 = [pydy_15; pydy_19; pydy_22; -input_4(1).*pydy_23 + ...
    input_4(2).*pydy_27 - input_4(8).*pydy_28 + pydy_10.*pydy_30 - ...
    pydy_22.*pydy_32 + pydy_31.*pydy_6; input_4(1).*pydy_15.*pydy_19 + ...
    input_4(1).*(-pydy_14.*pydy_16.*pydy_20 + pydy_16.*pydy_7 - ...
    pydy_18.*pydy_7 + pydy_20.*pydy_33 + pydy_25.*pydy_34) + ...
    input_4(2).*pydy_19.*pydy_22 - input_4(2).*(-pydy_16.*pydy_21 + ...
    pydy_16.*pydy_35 - pydy_18.*pydy_35 + pydy_25.*pydy_28 + ...
    pydy_33.*pydy_6) - input_4(8).*pydy_17 + pydy_36; -input_4(1).*pydy_27 ...
    - input_4(2).*pydy_23 - input_4(8).*pydy_34 + pydy_15.*pydy_32 + ...
    pydy_20.*pydy_31 + pydy_30.*pydy_6; pydy_47; pydy_54; pydy_59; ...
    input_4(4).*pydy_78 - input_4(5).*pydy_63 - input_4(6).*pydy_63 + ...
    input_4(7).*pydy_77 - input_4(8).*pydy_46 + pydy_36.*pydy_38 + ...
    pydy_38.*pydy_71 - pydy_40.*pydy_60 + pydy_40.*pydy_61 + ...
    pydy_40.*pydy_62 - pydy_59.*pydy_64 - pydy_59.*pydy_65; ...
    input_4(4).*pydy_80 - input_4(4).*(pydy_11.*(pydy_0.*pydy_51 - ...
    pydy_48.*pydy_69 + pydy_50.*pydy_55 + pydy_50.*pydy_56 + ...
    pydy_52.*pydy_55) + pydy_48.*pydy_55 - pydy_49.*pydy_55 + pydy_82) + ...
    input_4(5).*pydy_79 + input_4(5).*pydy_81 + input_4(6).*pydy_79 + ...
    input_4(6).*pydy_81 + input_4(7).*pydy_80 - ...
    input_4(7).*(pydy_19.*pydy_73 - pydy_48.*pydy_57 - pydy_67.*pydy_74 + ...
    pydy_72.*(pydy_3.*pydy_50 + pydy_51) + pydy_82) - input_4(8).*pydy_53 + ...
    pydy_36.*pydy_40 + pydy_38.*pydy_60 - pydy_38.*pydy_61 - ...
    pydy_38.*pydy_62 + pydy_40.*pydy_71; -input_4(3).*pydy_27 - ...
    input_4(4).*pydy_63 - input_4(5).*pydy_78 - input_4(6).*pydy_77 - ...
    input_4(7).*pydy_63 - input_4(8).*pydy_69 + pydy_3.*pydy_31 + ...
    pydy_30.*pydy_5 + pydy_47.*pydy_64 + pydy_47.*pydy_65 + ...
    pydy_58.*pydy_68];

end
