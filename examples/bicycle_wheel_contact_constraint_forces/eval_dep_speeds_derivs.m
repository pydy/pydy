function [output_1, output_2] = eval_dep_speeds_derivs(input_1, input_2, input_3, input_4)
% function [output_1, output_2] = eval_dep_speeds_derivs(input_1, input_2, input_3, input_4)
%
% input_1 : [q4(t), q5(t), q7(t)]
% input_2 : [u3(t), u4(t), u5(t), u6(t), u7(t), u8(t)]
% input_3 : [u4p(t), u6p(t), u7p(t)]
% input_4 : [d1, d2, d3, rf, rr]

    pydy_0 = input_1(1);
    pydy_1 = sin(pydy_0);
    pydy_2 = pydy_1.*input_4(5);
    pydy_3 = input_1(2);
    pydy_4 = sin(pydy_3);
    pydy_5 = input_4(1).*pydy_4;
    pydy_6 = pydy_1.*pydy_5;
    pydy_7 = input_1(3);
    pydy_8 = cos(pydy_7);
    pydy_9 = pydy_1.*pydy_8;
    pydy_10 = sin(pydy_7);
    pydy_11 = cos(pydy_0);
    pydy_12 = pydy_10.*pydy_11;
    pydy_13 = pydy_12.*pydy_4;
    pydy_14 = pydy_13 + pydy_9;
    pydy_15 = input_4(3).*pydy_14;
    pydy_16 = cos(pydy_3);
    pydy_17 = pydy_16.^2;
    pydy_18 = pydy_1.*pydy_10;
    pydy_19 = pydy_11.*pydy_8;
    pydy_20 = pydy_18.*pydy_4;
    pydy_21 = pydy_19 - pydy_20;
    pydy_22 = pydy_21.*pydy_4;
    pydy_23 = pydy_17.*pydy_18 - pydy_22;
    pydy_24 = pydy_4.*pydy_9;
    pydy_25 = pydy_12 + pydy_24;
    pydy_26 = pydy_10.*pydy_16;
    pydy_27 = pydy_16.*pydy_8;
    pydy_28 = pydy_21.*pydy_27 + pydy_25.*pydy_26;
    pydy_29 = pydy_23.^2 + pydy_28.^2;
    pydy_30 = input_4(4)./sqrt(pydy_29);
    pydy_31 = pydy_14.*pydy_30;
    pydy_32 = pydy_23.*pydy_31;
    pydy_33 = -pydy_15 - pydy_32;
    pydy_34 = input_4(2).*pydy_14 + pydy_28.*pydy_31;
    pydy_35 = pydy_19.*pydy_4;
    pydy_36 = pydy_18 - pydy_35;
    pydy_37 = input_4(2).*pydy_36;
    pydy_38 = pydy_11.*pydy_16;
    pydy_39 = input_4(3).*pydy_38;
    pydy_40 = pydy_23.*pydy_30;
    pydy_41 = pydy_38.*pydy_40;
    pydy_42 = pydy_28.*pydy_30;
    pydy_43 = pydy_36.*pydy_42;
    pydy_44 = -pydy_41 + pydy_43;
    pydy_45 = input_4(3).*pydy_8;
    pydy_46 = pydy_40.*pydy_8;
    pydy_47 = -pydy_45 - pydy_46;
    pydy_48 = input_4(2).*pydy_10;
    pydy_49 = pydy_10.*pydy_42;
    pydy_50 = input_4(2).*pydy_8 + pydy_42.*pydy_8;
    pydy_51 = input_4(1).*pydy_16;
    pydy_52 = pydy_11.^2.*pydy_51;
    pydy_53 = pydy_1.*pydy_16;
    pydy_54 = -pydy_37 + pydy_39 + pydy_41 - pydy_43;
    pydy_55 = pydy_1.*pydy_51;
    pydy_56 = -pydy_48 - pydy_49;
    pydy_57 = input_4(1).*pydy_38;
    pydy_58 = input_2(2);
    pydy_59 = input_2(1);
    pydy_60 = pydy_11.*pydy_59;
    pydy_61 = pydy_58.*pydy_60;
    pydy_62 = input_2(3);
    pydy_63 = pydy_1.*pydy_59;
    pydy_64 = pydy_62 + pydy_63;
    pydy_65 = pydy_62.*pydy_64;
    pydy_66 = pydy_62.*pydy_8;
    pydy_67 = pydy_26.*pydy_58;
    pydy_68 = pydy_14.*pydy_59 + pydy_66 - pydy_67;
    pydy_69 = input_4(3).*pydy_68;
    pydy_70 = pydy_68 + input_2(6);
    pydy_71 = pydy_40.*pydy_70;
    pydy_72 = -pydy_69 - pydy_71;
    pydy_73 = pydy_16.*pydy_62;
    pydy_74 = input_4(2).*pydy_68 + pydy_42.*pydy_70;
    pydy_75 = input_2(5);
    pydy_76 = pydy_26.*pydy_75;
    pydy_77 = pydy_4.*pydy_58;
    pydy_78 = pydy_38.*pydy_59 + pydy_77;
    pydy_79 = pydy_75 + pydy_78;
    pydy_80 = input_4(3).*pydy_79;
    pydy_81 = pydy_10.*pydy_62;
    pydy_82 = pydy_27.*pydy_58;
    pydy_83 = pydy_36.*pydy_59 + pydy_81 + pydy_82;
    pydy_84 = input_4(2).*pydy_83;
    pydy_85 = pydy_40.*pydy_79;
    pydy_86 = pydy_42.*pydy_83;
    pydy_87 = -pydy_80 + pydy_84 - pydy_85 + pydy_86;
    pydy_88 = pydy_4.*pydy_81;
    pydy_89 = pydy_27.*pydy_75;
    pydy_90 = input_3(1);
    pydy_91 = pydy_12.*pydy_73 - pydy_18.*pydy_75 - pydy_18.*pydy_77 + ...
    pydy_19.*pydy_58 + pydy_35.*pydy_75;
    pydy_92 = -pydy_26.*pydy_90 + pydy_59.*pydy_91 - pydy_75.*pydy_81 - ...
    pydy_75.*pydy_82 + pydy_77.*pydy_81;
    pydy_93 = input_4(3).*pydy_92;
    pydy_94 = pydy_40.*pydy_92;
    pydy_95 = pydy_12.*pydy_58;
    pydy_96 = pydy_17.*pydy_95;
    pydy_97 = pydy_75.*pydy_9;
    pydy_98 = pydy_17.*pydy_97;
    pydy_99 = pydy_20.*pydy_73;
    pydy_100 = pydy_73.*(-pydy_19 + pydy_20);
    pydy_101 = pydy_58.*pydy_9;
    pydy_102 = pydy_12.*pydy_75;
    pydy_103 = pydy_12.*pydy_77;
    pydy_104 = pydy_18.*pydy_73;
    pydy_105 = pydy_24.*pydy_75;
    pydy_106 = pydy_4.*(pydy_101 + pydy_102 + pydy_103 + pydy_104 + ...
    pydy_105);
    pydy_107 = pydy_100 + pydy_106 + pydy_96 + pydy_98 - 2*pydy_99;
    pydy_108 = pydy_30.*pydy_70;
    pydy_109 = pydy_107.*pydy_108;
    pydy_110 = pydy_25.*pydy_88;
    pydy_111 = pydy_25.*pydy_89;
    pydy_112 = pydy_22.*pydy_66;
    pydy_113 = pydy_21.*pydy_76;
    pydy_114 = -pydy_18.*pydy_58 + pydy_19.*pydy_75 + pydy_19.*pydy_77 - ...
    pydy_20.*pydy_75 + pydy_73.*pydy_9;
    pydy_115 = pydy_114.*pydy_26;
    pydy_116 = -pydy_101 - pydy_102 - pydy_103 - pydy_104 - pydy_105;
    pydy_117 = pydy_116.*pydy_27;
    pydy_118 = input_4(4).*(-pydy_23.*(2*pydy_100 + 2*pydy_106 + 2*pydy_96 ...
    + 2*pydy_98 - 4*pydy_99)/2 - pydy_28.*(-2*pydy_110 + 2*pydy_111 - ...
    2*pydy_112 - 2*pydy_113 + 2*pydy_115 + 2*pydy_117)/2)./pydy_29.^(3/2);
    pydy_119 = pydy_118.*pydy_70;
    pydy_120 = pydy_119.*pydy_23;
    pydy_121 = -pydy_109 - pydy_120 - pydy_93 - pydy_94;
    pydy_122 = -pydy_110 + pydy_111 - pydy_112 - pydy_113 + pydy_115 + ...
    pydy_117;
    pydy_123 = input_4(2).*pydy_92 + pydy_108.*pydy_122 + pydy_119.*pydy_28 ...
    + pydy_42.*pydy_92;
    pydy_124 = pydy_58.*pydy_63;
    pydy_125 = pydy_4.*pydy_62;
    pydy_126 = -pydy_124.*pydy_16 - pydy_125.*pydy_60 + pydy_4.*pydy_90 + ...
    pydy_58.*pydy_73;
    pydy_127 = pydy_126 + input_3(3);
    pydy_128 = input_4(3).*pydy_127;
    pydy_129 = pydy_13.*pydy_75 - pydy_19.*pydy_73 + pydy_77.*pydy_9 + ...
    pydy_95 + pydy_97;
    pydy_130 = pydy_129.*pydy_59 + pydy_27.*pydy_90 + pydy_66.*pydy_75 - ...
    pydy_66.*pydy_77 - pydy_67.*pydy_75;
    pydy_131 = input_4(2).*pydy_130;
    pydy_132 = pydy_127.*pydy_40;
    pydy_133 = pydy_107.*pydy_30.*pydy_79;
    pydy_134 = pydy_130.*pydy_42;
    pydy_135 = pydy_122.*pydy_30.*pydy_83;
    pydy_136 = pydy_118.*pydy_23.*pydy_79;
    pydy_137 = pydy_118.*pydy_28.*pydy_83;
    pydy_138 = pydy_11.*input_4(5);
    pydy_139 = pydy_58.^2;
    pydy_140 = pydy_58.*pydy_64;
    pydy_141 = input_4(1).*pydy_58.*pydy_78;
    pydy_142 = input_4(1).*pydy_126;
    pydy_143 = pydy_69 + pydy_71;
    pydy_144 = pydy_80 - pydy_84 + pydy_85 - pydy_86;
    pydy_145 = pydy_128 - pydy_131 + pydy_132 + pydy_133 - pydy_134 - ...
    pydy_135 + pydy_136 - pydy_137;

    output_1 = [-pydy_2 + pydy_26.*(pydy_37 - pydy_39 + pydy_44) + ...
    pydy_27.*pydy_34 + pydy_33.*pydy_4 - pydy_6 pydy_26.*(pydy_48 + ...
    pydy_49) + pydy_27.*pydy_50 + pydy_4.*pydy_47 - pydy_5 - input_4(5) ...
    pydy_27.*pydy_42 - pydy_4.*pydy_40; pydy_1.^2.*pydy_51 + ...
    pydy_21.*pydy_54 + pydy_25.*pydy_34 + pydy_52 + pydy_53.*(pydy_15 + ...
    pydy_32) pydy_21.*pydy_56 + pydy_25.*pydy_50 + pydy_53.*(pydy_45 + ...
    pydy_46) + pydy_55 pydy_25.*pydy_42 + pydy_40.*pydy_53; ...
    pydy_14.*pydy_54 + pydy_33.*pydy_38 + pydy_34.*pydy_36 pydy_14.*pydy_56 ...
    + pydy_36.*pydy_50 + pydy_38.*pydy_47 - pydy_57 pydy_44];

    output_2 = [-pydy_121.*pydy_4 - pydy_123.*pydy_27 - pydy_26.*(-pydy_128 ...
    + pydy_131 - pydy_132 - pydy_133 + pydy_134 + pydy_135 - pydy_136 + ...
    pydy_137) + pydy_4.*pydy_66.*pydy_74 + pydy_5.*pydy_61 + ...
    pydy_51.*pydy_65 - pydy_72.*pydy_73 + pydy_74.*pydy_76 + ...
    pydy_87.*pydy_88 - pydy_87.*pydy_89 + input_4(5).*(pydy_61 + ...
    input_3(2)); pydy_1.*pydy_125.*pydy_143 + pydy_1.*pydy_141 - ...
    pydy_11.*pydy_142 - pydy_114.*pydy_74 - pydy_116.*pydy_144 - ...
    pydy_123.*pydy_25 - pydy_124.*pydy_57 - pydy_138.*pydy_90 + ...
    pydy_139.*pydy_2 - pydy_140.*pydy_57 - pydy_143.*pydy_38.*pydy_58 - ...
    pydy_145.*pydy_21 - pydy_53.*(pydy_109 + pydy_120 + pydy_93 + pydy_94) ...
    + pydy_6.*pydy_65; -pydy_1.*pydy_142 + pydy_11.*pydy_125.*pydy_72 - ...
    pydy_11.*pydy_141 - pydy_11.*pydy_5.*pydy_65 - pydy_121.*pydy_38 - ...
    pydy_123.*pydy_36 - pydy_129.*pydy_74 - pydy_138.*pydy_139 - ...
    pydy_14.*pydy_145 - pydy_140.*pydy_55 - pydy_144.*pydy_91 - ...
    pydy_2.*pydy_90 + pydy_52.*pydy_58.*pydy_59 + ...
    pydy_53.*pydy_58.*pydy_72];

end
