function [legenda] = legendaFull
%LEGENDAFULL produces the legend for the different variables of study and
%writes them down in the variable legenda
%   Detailed explanation goes here

% legenda metabolites           %TO FILL TEMPLEG
legenda.metabolites = {'ACE, x_{1}, [mM]';
    'BPG, x_{2}, [mM]';
    'F16BP, x_{3}, [mM]';
    'F6P, x_{4}, [mM]';
    'G6P, x_{5}, [mM]';
    'GLCi, x_{6}, [mM]';
    'NAD, x_{7}, [mM]';
    'NADH, x_{8}, [mM]';
    'ATP, x_{9}, [mM]';
    'P2G, x_{10}, [mM]';
    'P3G, x_{11}, [mM]';
    'PEP, x_{12}, [mM]';
    'PYR, x_{13}, [mM]';
    'GLYCERAL3P, x_{14}, [mM]';
    'ADP, x_{15}, [mM]';
    'AMP, x_{16}, [mM]';
    'DHAP, x_{17}, [mM]';
    'GLYC3P, x_{18}, [mM]';
    'GLYCEROL, x_{19}, [mM]';
    'ETOH, x_{20}, [mM]';
    'G1P, x_{21}, [mM]';
    'UTP, x_{22}, [mM]';
    'UDP, x_{23}, [mM]';
    'UDP_{GLC}, x_{24}, [mM]';
    'TRE_{CYT}, x_{25}, [mM]';
    'T6P, x_{26}, [mM]';
    'PI, x_{27}, [mM]';
    'IMP, x_{28}, [mM]';
    'INO, x_{29}, [mM]';
    'HYP, x_{30}, [mM]';
    'ETOH_{EC}, x_{31}, [mM]';
    'GLYC_{EC}, x_{32}, [mM]';
    'FRC_{IC}, x_{33}, [mM]';
    'FRC_{EC}, x_{34}, [mM]';
    'SUC_{EC}, x_{35}, [mM]';
    'GLC_{EC}, x_{36}, [mM]';
    'TRE_{EC}, x_{37}, [mM]';
    'TRE_{vac}, x_{38}, [mM]'};

% legena reaction rates
legenda.fluxes = {'v_{GLT}, r_{1}, [mM s^{-1}]';
    'v_{GLK}, r_{2}, [mM s^{-1}]';
    'v_{PGI}, r_{3}, [mM s^{-1}]';
    'v_{PFK}, r_{4}, [mM s^{-1}]';
    'V_{ALD}, r_{5}, [mM s^{-1}]';
    'V_{G3PDH}, r_{6}, [mM s^{-1}]';
    'V_{TPI1}, r_{7}, [mM s^{-1}]';
    'V_{TDH1}, r_{8}, [mM s^{-1}]';
    'V_{PGK}, r_{9}, [mM s^{-1}]';
    'V_{PGM}, r_{10}, [mM s^{-1}]';
    'V_{ENO}, r_{11}, [mM s^{-1}]';
    'V_{PYK}, r_{12}, [mM s^{-1}]';
    'V_{PDC}, r_{13}, [mM s^{-1}]';
    'V_{ADK1}, r_{14}, [mM s^{-1}]';
    'V_{HOR2}, r_{15}, [mM s^{-1}]';
    'V_{RHR2}, r_{16}, [mM s^{-1}]';
    'V_{PGM1}, r_{17}, [mM s^{-1}]';
    'V_{UGP}, r_{18}, [mM s^{-1}]';
    'V_{TPS2}, r_{19}, [mM s^{-1}]';
    'V_{NTH1}, r_{20}, [mM s^{-1}]';
    'V_{TPS1}, r_{21}, [mM s^{-1}]';
    'V_{ACE}, r_{22}, [mM s^{-1}]';
    'V_{ETOHtransport}, r_{23}, [mM s^{-1}]';
    'V_{GLYCEROLtransport}, r_{24}, [mM s^{-1}]';
    'V_{PDH}, r_{25}, [mM s^{-1}]';
    'V_{mitoNADH}, r_{26}, [mM s^{-1}]';
    'V_{ATPase}, r_{27}, [mM s^{-1}]';
    'V_{mito}, r_{28}, [mM s^{-1}]';
    'V_{Amd1}, r_{29}, [mM s^{-1}]';
    'V_{Ade13_v_Ade12}, r_{30}, [mM s^{-1}]';
    'V_{Isn1}, r_{31}, [mM s^{-1}]';
    'V_{Pnp1}, r_{32}, [mM s^{-1}]';
    'V_{Hpt1}, r_{33}, [mM s^{-1}]';
    'V_{sinkG6P}, r_{34}, [mM s^{-1}]';
    'V_{sinkF6P}, r_{35}, [mM s^{-1}]';
    'V_{sinkGAP}, r_{36}, [mM s^{-1}]';
    'V_{sinkP3G}, r_{37}, [mM s^{-1}]';
    'V_{sinkPEP}, r_{38}, [mM s^{-1}]';
    'V_{sinkPYR}, r_{39}, [mM s^{-1}]';
    'V_{sinkACE}, r_{40}, [mM s^{-1}]';
    'V_{ADH}, r_{41}, [mM s^{-1}]';
    'V_{vacuolePi}, r_{42}, [mM s^{-1}]';
    'V_{FRT}, r_{43}, [mM s^{-1}]';
    'V_{FRK}, r_{44}, [mM s^{-1}]';
    'V_{ATH1ec}, r_{45}, [mM s^{-1}]';
    'V_{ATH1vac}, r_{46}, [mM s^{-1}]';
    'V_{AGT1}, r_{47}, [mM s^{-1}]';
    'V_{vacuoleT}, r_{48}, [mM s^{-1}]'};

% legenda parameters
legenda.parameters = cell(1,156);

%ADH
legenda.parameters(1) = {'px1, Keq, mM'};
legenda.parameters(2) = {'px2, KiACE, mM'};
legenda.parameters(3) = {'px3, KiETOH, mM'};
legenda.parameters(4) = {'px4, KiNAD, mM'};
legenda.parameters(5) = {'px5, KiNADH, mM'};
legenda.parameters(6) = {'px6, KmACE, mM'};
legenda.parameters(7) = {'px7, KmETOH, mM'};
legenda.parameters(8) = {'px8, KmNAD, mM'};
legenda.parameters(9) = {'px9, KmNADH, mM'};
legenda.parameters(10) = {'px10, Vm, s^{-1}'};

%FBA
legenda.parameters(11) = {'px11, kdhap, mM'};
legenda.parameters(12) = {'px12, keq, []'};
legenda.parameters(13) = {'px13, kf16bp, mM'};
legenda.parameters(14) = {'px14, kgap, mM'};
legenda.parameters(15) = {'px15, kcat, mM s^{-1}'};

%ENO
legenda.parameters(16) = {'px16, k2pg, mM'};
legenda.parameters(17) = {'px17, keq, []'};
legenda.parameters(18) = {'px18, kpep, mM'};
legenda.parameters(19) = {'px19, kcat, mM s^{-1}'};

%TDH1/GAPDH
legenda.parameters(21) = {'px21, keq, []'};
legenda.parameters(22) = {'px22, kgap, mM'};
legenda.parameters(23) = {'px23, kbpg, mM'};
legenda.parameters(24) = {'px24, knad, mM'};
legenda.parameters(25) = {'px25, knadh, mM'};
legenda.parameters(26) = {'px26, kpi, mM'};
legenda.parameters(27) = {'px27, kcat, mM s^{-1}'};

%GLK
legenda.parameters(28) = {'px_{28}, k_{ADP}, mM'};
legenda.parameters(29) = {'px_{29}, k_{ATP}, mM'};
legenda.parameters(30) = {'px_{30}, k_{eq}, []'};
legenda.parameters(31) = {'px_{31}, k_{G6P}, mM'};
legenda.parameters(32) = {'px_{32}, k_{GLC}, mM'};
legenda.parameters(33) = {'px_{33}, k_{T6P}, mM'};
legenda.parameters(34) = {'px_{34}, k_{cat}, s^{-1}'};

%GLT
legenda.parameters(35) = {'px_{35}, k_{GLCoGLCi}, mM'};
legenda.parameters(36) = {'px_{36}, r_{m}, mM'};

%PDC
legenda.parameters(39) = {'px39, kpi, mM'};
legenda.parameters(40) = {'px40, kpyr, mM'};
legenda.parameters(41) = {'px41, hill, []'};
legenda.parameters(42) = {'px42, kcat, mM s^{-1}'};

%PFK
legenda.parameters(43) = {'px_{43}, c_{iAMP}, []'};
legenda.parameters(44) = {'px_{44}, c_{ATP}, []'};
legenda.parameters(45) = {'px_{45}, c_{F16BP}, []'};
legenda.parameters(46) = {'px_{46}, c_{F26BP}, [],'};
legenda.parameters(47) = {'px_{47}, c_{ATP}, []'};
legenda.parameters(48) = {'px_{48}, k_{AMP}, mM'};
legenda.parameters(49) = {'px_{49}, k_{ATP}, mM'};
legenda.parameters(50) = {'px_{50}, k_{F16BP}, mM'};
legenda.parameters(51) = {'px_{51}, k_{F26BP}, mM'};
legenda.parameters(52) = {'px_{52}, k_{F6P}, mM'};
legenda.parameters(53) = {'px_{53}, k_{iATP}, mM'};
legenda.parameters(54) = {'px_{54}, L, []'};
legenda.parameters(55) = {'px_{55}, gR, []'};
legenda.parameters(56) = {'px_{56}, k_{cat}, s^{-1}'};
legenda.parameters(87) = {'px_{87}, F26BP, mM'};

% PGI
legenda.parameters(57) = {'px_{57}, k_{eq}, []'};
legenda.parameters(58) = {'px_{58}, k_{F6P}, mM'};
legenda.parameters(59) = {'px_{59}, k_{G6P}, mM'};
legenda.parameters(60) = {'px_{60}, k_{cat}, s^{-1}'};

%TPI
legenda.parameters(57) = {'px57, Keq, []'};
legenda.parameters(58) = {'px58, Kf6p, mM'};
legenda.parameters(59) = {'px59, Kg6p, mM'};
legenda.parameters(60) = {'px60, kcat, mM s^{-1}'};

%PGK
legenda.parameters(61) = {'px61, keq, []'};
legenda.parameters(62) = {'px62, kadp, mM'};
legenda.parameters(63) = {'px63, katp, mM'};
legenda.parameters(64) = {'px64, kbpg, mM'};
legenda.parameters(65) = {'px65, kp3g, mM'};
legenda.parameters(66) = {'px66, rmax, mM s^{-1}'};

%PGM (GPM1)
legenda.parameters(67) = {'px67, k2pg, mM'};
legenda.parameters(68) = {'px68, k3pg, mM'};
legenda.parameters(69) = {'px69, keq, []'};
legenda.parameters(70) = {'px70, kcat, mM s^{-1}'};

%PYK
legenda.parameters(71) = {'px71, kadp, mM'};
legenda.parameters(72) = {'px72, katp, mM'};
legenda.parameters(73) = {'px73, kf16bp, mM'};
legenda.parameters(74) = {'px74, kpep, mM'};
legenda.parameters(75) = {'px75, L, []'};
legenda.parameters(76) = {'px76, hill, []'};
legenda.parameters(77) = {'px77, kcat, mM s^{-1}'};

%TPI
legenda.parameters(79) = {'px79, kdhap, mM'};
legenda.parameters(80) = {'px80, keq, []'};
legenda.parameters(81) = {'px81, kgap, mM'};
legenda.parameters(82) = {'px82, kcat, mM s^{-1}'};

%PGM1
legenda.parameters(83) = {'px83, keq^{PGM1}, '};
legenda.parameters(84) = {'px84, kg1p^{PGM1}, '};
legenda.parameters(85) = {'px85, kg6p^{PGM1}, '};
legenda.parameters(86) = {'px86, kcat^{PGM1}, '};

%GPD (G3PDH)
legenda.parameters(96)  = {'px96, kadp, mM'};
legenda.parameters(97)  = {'px97, katp, mM'};
legenda.parameters(98)  = {'px98, kdhap, mM'};
legenda.parameters(99)  = {'px99, keq, []'};
legenda.parameters(100) = {'px100, kf16bp, mM'};
legenda.parameters(101) = {'px101, kglyc3p, mM'};
legenda.parameters(102) = {'px102, knad, mM'};
legenda.parameters(103) = {'px103, knadh, mM'};
legenda.parameters(104) = {'px104, kcat, mM s^{-1}'};

%HOR2 (GlycT)
legenda.parameters(105) = {'px105, kglyc3p, '};
legenda.parameters(106) = {'px106, Kpi, '};
legenda.parameters(107) = {'px107, kcat, '};
legenda.parameters(108) = {'px108, KGlycerolTransport, '};

% % % % % Trehalose cycle % % % % %

%TPS1
legenda.parameters(124) = {'px124, k_{m,G6P}^{TPS1}, mM'};
legenda.parameters(125) = {'px125, k_{m,udpglc}^{TPS1}, mM'};
legenda.parameters(126) = {'px126, k_{cat}^{TPS1}, s^{-1}'};
legenda.parameters(127) = {'px127, k_{m,pi}^{TPS1}, '};
legenda.parameters(128) = {'px128, k_{m,F6P}^{TPS1}, '};
    
%TPS2
legenda.parameters(119) = {'px119, k^{TPS2}_{m,t6p}, mM'};
legenda.parameters(120) = {'px120, k^{TPS2}_{cat}, s^{-1}'};
legenda.parameters(121) = {'px121, k^{TPS2}_{m,pi}, '};
    
%NTH1
legenda.parameters(122) = {'px122, k^{NTH1}_{m,tre}, mM'};
legenda.parameters(123) = {'px123, k^{NTH1}_{cat}, s^{-1}'};

%UGP
legenda.parameters(144) = {'px144, k^{UGP}_{m,utp}, mM'};
legenda.parameters(145) = {'px145, k^{UGP}_{i,utp}, mM'};
legenda.parameters(146) = {'px146, k^{UGP}_{m,g1p}, mM'};
legenda.parameters(147) = {'px147, k^{UGP}_{i,UDPgluc}, mM'};
legenda.parameters(148) = {'px148, k^{UGP}_{cat}, s^{-1}'};

%ATH1
legenda.parameters(149) = {'px149, k^{ATH1}_{m,tre}, mM'};
legenda.parameters(150) = {'px150, kcat_{ATH1}, mM/s^{-1}'};

%AGT1
legenda.parameters(152) = {'px152, kcat_{AGT1}, s^{-1}'};
legenda.parameters(153) = {'px153, k^{AGT1}_{m,TREic}, mM'};
legenda.parameters(154) = {'px154, k^{AGT1}_{eq}'};
legenda.parameters(151) = {'px151, k^{AGT1}_{m,TREec}, mM'};
legenda.parameters(158) = {'px158, K^{AGT1}_{i, UDPG}'};

%Transport vacuoleT
legenda.parameters(155) = {'px155, vmax_{vacuoleT}, mM/s^{-1}'};
legenda.parameters(156) = {'px156, k^{vacuoleT}_{m,tre}, mM'};
legenda.parameters(157) = {'px157, K^{vacuoleT}_{eq}'};

% % % % % Cofactors % % % % %

%mitoNADH
legenda.parameters(92) = {'px92, v^{mitoNADH}_{max}, mM'};
legenda.parameters(93) = {'px93, k^{mitoNADH}_{m,nadh}, mM'};

%vacImport
legenda.parameters(37) = {'px37, k^{vacPi}_{m}, mM^{-1} s^{-1}'};
legenda.parameters(38) = {'px38, Pi^{vacPi}_{SS}, mM'};

%mito: mitochondrial ATP synthesis
legenda.parameters(109) = {'px109, v^{mito}_{max}, '};
legenda.parameters(110) = {'px110, k^{mito}_{m,adp}, '};
legenda.parameters(111) = {'px111, k^{mito}_{m,pi}, '};

%ATPase: cytoplasmic ATPase
legenda.parameters(129) = {'px129, k^{ATPase}_{k}, '};

%ADK: ADK1 consensus
legenda.parameters(130) = {'px130, k^{ADK1}_{keq}, '};
legenda.parameters(131) = {'px131, k^{ADK1}_{k}, '};

% Fructose transport
legenda.parameters(137) = {'px137, KmFRTFRCi, mM'};
legenda.parameters(138) = {'px138, KmFRTFRCo, mM'};
legenda.parameters(139) = {'px139, vMax FRT, mM/s'};

%HXK Fructose
legenda.parameters(140) = {'px140, HXK1 Keq_{frc}'};
legenda.parameters(141) = {'px141, HXK1 Kf6p'};
legenda.parameters(142) = {'px142, HXK1 Kfrc'};
legenda.parameters(143) = {'px143,HXK1 kcatFrc'};

end

