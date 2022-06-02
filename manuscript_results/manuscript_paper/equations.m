


% latex(...
%     dACE/dt == 1...
%     )

syms dACE dt v_PDC v_ADH v_sinkACE
latex(...
    dACE/dt == + v_PDC - v_ADH - v_sinkACE)

syms dACE_L dt v_PDC v_ADH v_sinkACE PYR_L ACE_L
latex(...
    dACE_L/dt == + PYR_L * v_PDC - ACE_L * v_ADH - ACE_L * v_sinkACE)

%% v_GLT

syms v_GLT V_m GLC_ec2 GLC_i K_eq K_m K_i
latex(...
    v_GLT == V_m*(GLC_ec2-GLC_i/K_eq)/...
    (K_m*(1+GLC_ec2/K_m+GLC_i / K_m + ...
    K_i*GLC_ec2*GLC_i/(K_m^2)))...
    )
syms GLC_ec2 GLC_ec GLC_ec_min
latex(...
    GLC_ec2 == GLC_ec - GLC_ec_min...
    )
syms GLC_ec GLC_ec_min
latex(...
    GLC_ec > GLC_ec_min...
    )


%% regularization
syms error_estimation error_data error_parameters
latex(error_estimation == error_data + error_parameters)
syms error_data data_simulated data_experimental
latex(error_data == data_simulated - data_experimental)
syms error_parameters lambda parameters_estimated parameters_reference
latex(error_parameters == lambda * (parameters_estimated - parameters_reference))


%% Y3M2 added reactions

syms v_ATH K_cat ATH TRE K_M_TRE
latex(...
    v_ATH == (K_cat * ATH)* ( (TRE/K_M_TRE)/ ( 1+(TRE/K_M_TRE) ) )...
    )

syms v_AGT1 K_cat AGT1 K_M_TRE TRE_cyt TRE_ec K_eq UDP_GLC K_i_UDPglc
latex(...
    v_AGT1 == (K_cat * AGT1)* (1/ K_M_TRE) * ( TRE_cyt - TRE_ec/K_eq )/ ( 1+TRE_cyt/K_M_TRE + TRE_ec/K_M_TRE + (UDP_GLC/K_i_UDPglc) ) ...
    )

syms v_vacuoleT V_max K_M_TRE TRE_cyt TRE_vac K_eq
latex(...
    v_vacuoleT == V_max * (1/K_M_TRE) * ( TRE_cyt - TRE_vac/K_eq) / ( 1 + TRE_cyt/ K_M_TRE + TRE_vac/ K_M_TRE)...
    )

%% Glycogen
syms v_glyc_synthesis v_glyc_synthesis_interpolated UDP_GLC
latex(...
    v_glyc_synthesis == v_glyc_synthesis_interpolated * UDP_GLC / (UDP_GLC + 0.0001)...
    )
syms v_glyc_degradation v_glyc_degradation_interpolated Glycogen
latex(...
    v_glyc_degradation == v_glyc_degradation_interpolated  * Glycogen / (Glycogen + 0.0001)...
    )


    
% % converting equations to latex format
% 
% % % sample
% % syms x phi
% % latex(x^2 + 1/x)
% % latex(x^2 == 1/x)
% % latex(x^2 == ...
% %     1/x)
% % syms v_GLT
% % latex(...
% %     v_GLT...
% %     )
% % syms k_m_GLT
% % latex(...
% %     k_m_GLT...
% %     )
% 
% %% v_GLT
% % v_GLT=p.GLT.VmGLT*(f.GLCo-GLCi/p.GLT.KeqGLT)/...
% %     (p.GLT.KmGLTGLCo*(1+f.GLCo/p.GLT.KmGLTGLCo+GLCi / p.GLT.KmGLTGLCi + ...
% %     0.91*f.GLCo*GLCi/(p.GLT.KmGLTGLCi*p.GLT.KmGLTGLCo)));
% syms v_GLT V_m GLCo GLCi K_eq K_m K_i
% latex(...
%     v_GLT == V_m*(GLCo-GLCi/K_eq)/...
%     (K_m*(1+GLCo/K_m+GLCi / K_m + ...
%     K_i*GLCo*GLCi/(K_m^2)))...
%     )
% 
% %% v_HXK
% % v_GLK=p.HXK_ExprsCor*(((p.HXK1_kcat*(f.HXK1+f.HXK2))/...
% %     (p.HXK1_Katp*p.HXK1_Kglc)*(ATP*GLCi-((ADP*G6P)/p.HXK1_Keq)))/...
% %     ((1+ATP/p.HXK1_Katp+ADP/p.HXK1_Kadp)*(1+GLCi/p.HXK1_Kglc+G6P/...
% %     p.HXK1_Kg6p+T6P/p.HXK1_Kt6p)));
% syms v_HXK HXK_cor V_m K_m_atp K_m_glc ATP GLCi ADP G6P K_eq K_m_adp K_m_g6p T6P K_i_t6p
% latex(v_HXK == ((V_m*HXK_cor/...
%     (K_m_atp*K_m_glc)*(ATP*GLCi-((ADP*G6P)/K_eq)))/...
%     ((1+ATP/K_m_atp+ADP/K_m_adp)*(1+GLCi/K_m_glc+G6P/...
%     K_m_g6p+T6P/K_i_t6p)))...
%     )
% 
% 
% %% PGI1 
% % v_PGI=p.PGI_ExprsCor*((((p.PGI1_kcat*f.PGI1)/p.PGI1_Kg6p)*(G6P-(F6P/p.PGI1_Keq)))/...
% %     (1+G6P/p.PGI1_Kg6p+1+F6P/p.PGI1_Kf6p5));
% syms v_PGI V_m PGI_cor K_m_g6p K_m_f6p K_eq G6P F6P
% latex(v_PGI == ((((V_m*PGI_cor)/...
%     K_m_g6p)*(G6P-(F6P/K_eq)))/...
%     (1+G6P/K_m_g6p+F6P/K_m_f6p))...
%     )
% 
% %% v_PFK
% % F26BP=p.PFK.F26BP;
% % PFK_nom=(p.PFK_kcat*f.PFK*p.PFK_gR*(F6P/p.PFK_Kf6p)*(ATP/p.PFK_Katp)*(1+(F6P/p.PFK_Kf6p)+(ATP/p.PFK_Katp)+p.PFK_gR*((F6P/p.PFK_Kf6p)*(ATP/p.PFK_Katp))));
% % PFK_denom=(1+F6P/p.PFK_Kf6p+ATP/p.PFK_Katp+(p.PFK_gR*(F6P/p.PFK_Kf6p)*(ATP/p.PFK_Katp))).^2+...
% %     p.PFK_L*...
% %     ((1+p.PFK_Ciatp*(ATP/p.PFK_Kiatp))/(1+ATP/p.PFK_Kiatp)).^2*...
% %     ((1+p.PFK_Camp*(AMP/p.PFK_Kamp))/(1+AMP/p.PFK_Kamp)).^2*...
% %     ((1+((p.PFK_Cf26bp*F26BP)/(p.PFK_Kf26bp))+((p.PFK_Cf16bp*F16BP)/(p.PFK_Kf16bp)))/(1+(F26BP/p.PFK_Kf26bp)+(F16BP/p.PFK_Kf16bp))).^2*...
% %     (1+p.PFK_Catp*(ATP/p.PFK_Katp)).^2;
% % v_PFK=p.PFK_ExprsCor*(PFK_nom/PFK_denom);
% syms v_PFK v_m PFK_cor g_R lam1 lam2 R L T
% latex(v_PFK == v_m * PFK_cor*...
%     (g_R * lam1 * lam2 * R)/...
%     (R^2 + L * T^2)...
%     )
% % %% 
% syms lam1 F6P K_R_F6P
% latex(lam1 == F6P/K_R_F6P)
% % %% 
% syms lam2 ATP K_R_ATP
% latex(lam2 == ATP/K_R_ATP)
% % %%
% syms R lam1 lam2 g_R
% latex(R == 1 + lam1 * lam2 + lam1 * lam2 * g_R)
% % %%
% syms T c_ATP lam2
% latex(T == 1 + c_ATP * lam2)
% % %%
% syms L L0 c_i_ATP c_i_AMP ATP AMP K_ATP K_AMP c_i_F26bp c_i_FBP F26bP FBP K_F26bP K_FBP
% latex(L == L0 * ...
%     ((1 + c_i_ATP * ATP / K_ATP)/(1 + ATP / K_ATP))^2 * ...
%     ((1 + c_i_AMP * AMP / K_AMP)/(1 + AMP / K_AMP))^2 * ...
%     ((1 + c_i_F26bp * F26bP / K_F26bP + c_i_FBP * FBP / K_FBP)/...
%     (1 + F26bP / K_F26bP + FBP / K_FBP))...
%     )
% 
% %% FBA1/ALD
% syms v_ALD FBA_cor V_m K_m_FBP FBP GAP DHAP K_eq K_m_gap K_m_dhap
% latex(v_ALD == ((V_m*FBA_cor/...
%     K_m_FBP*(FBP-(GAP*DHAP)/K_eq))/...
%     (1+FBP/K_m_FBP+(1+GAP/K_m_gap)*...
%     (1+DHAP/K_m_dhap)-1))...
%     )
% % v_ALD=p.FBA_ExprsCor*(((p.FBA1_kcat*f.FBA1)/p.FBA1_Kf16bp*(F16BP-(GLYCERAL3P*DHAP)/p.FBA1_Keq))/...
% %     (1+F16BP/p.FBA1_Kf16bp+(1+GLYCERAL3P/p.FBA1_Kglyceral3p)*(1+DHAP/p.FBA1_Kdhap)-1));
% 
% %% TPI1
% syms v_TPI TPI_cor V_m K_m_dhap DHAP GAP K_eq K_m_gap
% latex(v_TPI == (((V_m*TPI_cor)/K_m_dhap*(DHAP-GAP/K_eq))/...
%     (1+DHAP/K_m_dhap+1+GAP/K_m_gap-1))...
%     )
% % v_TPI1=(((p.TPI1_kcat*f.TPI1)/p.TPI1_Kdhap*(DHAP-GLYCERAL3P/p.TPI1_Keq))/...
% %     (1+DHAP/p.TPI1_Kdhap+1+GLYCERAL3P/p.TPI1_Kglyceral3p-1));
% 
% %% TDH1/GAPDH 
% syms v_GAPDH GAPDH_cor V_m K_m_gap K_m_nad K_m_pi GAP NAD PI BPG NADH K_eq K_m_bpg K_m_nadh
% latex(v_GAPDH == GAPDH_cor*(((V_m/...
%     (K_m_gap * K_m_nad * K_m_pi))*(GAP*NAD*...
%     PI-(BPG*NADH)/K_eq))/...
%     ((1+GAP/K_m_gap)*(1+NAD/K_m_nad)*...
%     (1+PI/K_m_pi)+(1+BPG/K_m_bpg)*...
%     (1+NADH/K_m_nadh)-1))...
%     )
% v_GAPDH=p.GAPDH_ExprsCor*((((p.TDH1_kcat*(f.TDH1+f.TDH2+f.TDH3))/(p.TDH1_Kglyceral3p*p.TDH1_Knad*p.TDH1_Kpi))*(GLYCERAL3P*NAD*PI-(BPG*NADH)/p.TDH1_Keq))/...
%     ((1+GLYCERAL3P/p.TDH1_Kglyceral3p)*(1+NAD/p.TDH1_Knad)*(1+PI/p.TDH1_Kpi)+(1+BPG/p.TDH1_Kglycerate13bp)*(1+NADH/p.TDH1_Knadh)-1));
% 
% %% PGK teusink
% syms v_PGK PGK_cor V_m K_eq BPG ADP ATP P3G K_m_ATP K_m_P3G K_m_ADP K_m_BPG
% latex(v_PGK == PGK_cor*V_m*((K_eq*BPG*ADP)-ATP*P3G)/...
%   (K_m_ATP*K_m_P3G*(1+ADP/K_m_ADP + ATP/K_m_ATP)*(1+BPG/K_m_BPG+P3G/K_m_P3G))...
%   )
% % v_PGK=p.PGK_ExprsCor*p.PGK.VmPGK*((p.PGK.KeqPGK*BPG*ADP)-ATP*P3G)/...
% %   (p.PGK.KmPGKATP*p.PGK.KmPGKP3G*(1+ADP/p.PGK.KmPGKADP + ATP/p.PGK.KmPGKATP)*(1+BPG/p.PGK.KmPGKBPG+P3G/p.PGK.KmPGKP3G));
%   
% %% GPM1 Consensus
% syms v_PGM PGM_cor V_m K_m_P3G K_m_P2G P3G P2G K_eq
% latex(v_PGM == PGM_cor*(((V_m/K_m_P3G)*(P3G-P2G/K_eq))/...
%     (1+P3G/K_m_P3G+1+P2G/K_m_P2G-1))...
%     )
% % v_PGM=p.PGM_ExprsCor*((((p.GPM1_kcat*(f.GPM1+f.GPM2+f.GPM3))/p.GPM1_K3pg)*(P3G-P2G/p.GPM1_Keq))/...
% %     (1+P3G/p.GPM1_K3pg+1+P2G/p.GPM1_K2pg-1));
% 
% %% ENO1
% syms v_ENO ENO_cor V_m K_m_PEP K_m_P2G PEP P2G K_eq
% latex(v_ENO == ENO_cor*(((V_m/K_m_P2G)*(P2G-PEP/K_eq))/...
%     (1+P2G/K_m_P2G+1+PEP/K_m_PEP-1))...
%     )
% % v_ENO=p.ENO_ExprsCor*((((p.ENO1_kcat*(f.ENO1+f.ENO2))/p.ENO1_K2pg)*(P2G-PEP/p.ENO1_Keq))/...
% %     (1+P2G/p.ENO1_K2pg+1+PEP/p.ENO1_Kpep-1));
% 
% %% PYK1 Consensus
% syms v_PYK PYK_cor V_m K_m_adp K_m_pep ADP PEP hill L K_a_FBP K_i_ATP FBP ATP
% latex(v_PYK == PYK_cor*(((V_m/(K_m_adp*K_m_pep)*ADP*PEP)/...
%     ((1+ADP/K_m_adp)*(1+PEP/K_m_pep)))*...
%     ((PEP/K_m_pep+1).^hill/(L*((ATP/K_i_ATP+1)/...
%     (FBP/K_a_FBP+1)).^hill+(PEP/K_m_pep+1).^hill)))...
%     )
% % v_PYK=p.PYK_ExprsCor*((((p.PYK1_kcat*(f.PYK1+f.PYK2))/(p.PYK1_Kadp*p.PYK1_Kpep)*ADP*PEP)/...
% %     ((1+ADP/p.PYK1_Kadp)*(1+PEP/p.PYK1_Kpep)))*...
% %     ((PEP/p.PYK1_Kpep+1).^p.PYK1_hill/(p.PYK1_L*((ATP/p.PYK1_Katp+1)/(F16BP/p.PYK1_Kf16bp+1)).^p.PYK1_hill+(PEP/p.PYK1_Kpep+1).^p.PYK1_hill)));
% 
% %% PDC
% syms v_PDC PDC_cor V_m PYR K_m_PYR hill 
% latex(v_PDC == PDC_cor*((V_m*(PYR/K_m_PYR)^hill)/...
%     (1+(PYR/K_m_PYR)^hill))...
%     )
% % v_PDC=p.PDC_ExprsCor*((p.PDC1_kcat*(f.PDC1)*(PYR/p.PDC1_Kpyr).^p.PDC1_hill)/...
% %     (1+(PYR/p.PDC1_Kpyr).^p.PDC1_hill));%+PI/p.PDC1_Kpi));
% 
% %% ADH 
% syms v_ADH ADH_cor V_m K_i_NAD K_i_NADH K_m_ETOH NADH ACE K_eq NAD ETOH K_m_NAD K_m_NADH K_m_ACE K_i_ETOH K_i_ACE
% latex(v_ADH == ADH_cor*(V_m/(K_i_NAD*K_m_ETOH)*...
%     (NADH*ACE/K_eq - NAD*ETOH)/...
%     (1+NAD/K_i_NAD+K_m_NAD*ETOH/(K_i_NAD*K_m_ETOH)...
%     +K_m_NADH*ACE/(K_i_NADH*K_m_ACE)...
%     +NADH/K_i_NADH+NAD*ETOH/(K_i_NAD*K_m_ETOH)...
%     +K_m_NADH*NAD*ACE/(K_i_NAD*K_i_NADH*K_m_ACE)...
%     +K_m_NAD*ETOH*NADH/(K_i_NAD*K_m_ETOH*K_i_NADH)...
%     +NADH*ACE/(K_i_NADH*K_m_ACE)+...
%     NAD*ETOH*ACE/(K_i_NAD*K_m_ETOH*K_i_ACE)...
%     +ETOH*NADH*ACE/(K_i_ETOH*K_i_NADH*K_m_ACE)))...
%     )
% % v_ADH=-p.ADH_ExprsCor*(p.ADH.VmADH/(p.ADH.KiADHNAD*p.ADH.KmADHETOH)*(NAD*ETOH-NADH*ACE/p.ADH.KeqADH)/...
% %     (1+NAD/p.ADH.KiADHNAD+p.ADH.KmADHNAD*ETOH/(p.ADH.KiADHNAD*p.ADH.KmADHETOH)+p.ADH.KmADHNADH*ACE/(p.ADH.KiADHNADH*p.ADH.KmADHACE)...
% %     +NADH/p.ADH.KiADHNADH+NAD*ETOH/(p.ADH.KiADHNAD*p.ADH.KmADHETOH)+p.ADH.KmADHNADH*NAD*ACE/(p.ADH.KiADHNAD*p.ADH.KiADHNADH*p.ADH.KmADHACE)...
% %     +p.ADH.KmADHNAD*ETOH*NADH/(p.ADH.KiADHNAD*p.ADH.KmADHETOH*p.ADH.KiADHNADH)+NADH*ACE/(p.ADH.KiADHNADH*p.ADH.KmADHACE)+...
% %     NAD*ETOH*ACE/(p.ADH.KiADHNAD*p.ADH.KmADHETOH*p.ADH.KiADHACE)+ETOH*NADH*ACE/(p.ADH.KiADHETOH*p.ADH.KiADHNADH*p.ADH.KmADHACE)));
% 
% %% GPD1
% syms v_GPD V_m K_m_DHAP K_m_NADH DHAP NADH G3P NAD K_eq FBP K_i_FBP ATP K_i_ATP ADP K_i_ADP K_m_G3P K_m_NAD
% latex(v_GPD == (((V_m/(K_m_DHAP*K_m_NADH))*(DHAP*NADH-(G3P*NAD)/K_eq))/...
%     ((1+FBP/K_i_FBP+ATP/K_i_ATP+ADP/K_i_ADP)*...
%     (1+DHAP/K_m_DHAP+G3P/K_m_G3P)*(1+NADH/K_m_NADH+NAD/K_m_NAD)))...
%     )
% % v_G3PDH=((((p.GPD1_kcat*f.GPD1)/(p.GPD1_Kdhap*p.GPD1_Knadh))*(DHAP*NADH-(GLYC3P*NAD)/p.GPD1_Keq))/...
% %     ((1+F16BP/p.GPD1_Kf16bp+ATP/p.GPD1_Katp+ADP/p.GPD1_Kadp)*(1+DHAP/p.GPD1_Kdhap+GLYC3P/p.GPD1_Kglyc3p)*(1+NADH/p.GPD1_Knadh+NAD/p.GPD1_Knad)));
% 
% %% HOR2 
% syms v_HOR2 V_m K_m_G3P G3P PI K_i_PI 
% latex(v_HOR2 == ((V_m/K_m_G3P*G3P)/...
%     ((1+PI/K_i_PI)*(1+G3P/K_m_G3P)))...
%     )
% % v_HOR2=(((p.HOR2_kcat*f.HOR2)/p.HOR2_Kglyc3p*GLYC3P)/...
% %     ((1+PI/p.HOR2_Kpi)*(1+GLYC3P/p.HOR2_Kglyc3p)));
% 
% %% PGM1
% syms v_PGM1 V_m K_m_G1P G1P G6P K_eq K_m_G6P
% latex(v_PGM1 == ((V_m/K_m_G1P*(G1P-G6P/K_eq))/...
%     (1+G1P/K_m_G1P+G6P/K_m_G6P))...
%     )
% % v_PGM1=(((p.PGM1_kcat*(f.PGM1+f.PGM2+f.PGM3))/p.PGM1_Kg1p*(G1P-G6P/p.PGM1_Keq))/...
% %     (1+G1P/p.PGM1_Kg1p+G6P/p.PGM1_Kg6p));
% 
% %% TPS1
% syms v_TPS1 F6P K_m_F6P V_m K_m_G6P K_m_UDPGLC G6P UDPGLC K_m_UDPGLC PI K_i_PI
% % 
% latex(v_TPS1 == (F6P/(F6P+K_m_F6P))*((V_m/...
%     (K_m_G6P*K_m_UDPGLC)*G6P*UDPGLC/...
%     ((1+G6P/K_m_G6P)*(1+UDPGLC/...
%     K_m_UDPGLC)*(1+PI/K_i_PI))))...
%     )
% % v_TPS1=(F6P/(F6P+p.TPS1_KmF6P))*(((p.TPS1_kcat*f.TPS1)/(p.TPS1_Kg6p*p.TPS1_Kudp_glc)*G6P*UDP_GLC/...
% %     ((1+G6P/p.TPS1_Kg6p)*(1+UDP_GLC/p.TPS1_Kudp_glc)*(1+PI/p.TPS1_Kpi))));
% 
% 
% %% TPS2
% syms v_TPS2 V_m T6P PI K_m_T6P K_i_PI
% latex(v_TPS2 == ((V_m*T6P*PI)/...
%     ((K_m_T6P*K_i_PI)+(K_m_T6P+T6P)*PI))...
%     )
% % v_TPS2=(((p.TPS2_kcat*f.TPS2)*T6P*PI)/...
% %     ((p.TPS2_Kt6p*p.TPS2_Kpi)+(p.TPS2_Kt6p+T6P)*PI));
% 
% %% NTH1
% syms v_NTH1 V_m K_m_TRE TRE
% latex(v_NTH1==((V_m/K_m_TRE*TRE)/...
%     (1+TRE/K_m_TRE))...
%     )
% % v_NTH1=(((p.NTH1_kcat*f.NTH1)/p.NTH1_Ktre*TRE)/...
% %     (1+TRE/p.NTH1_Ktre));
% 
% %% ethanol and glycerol transport
% syms V_ETOHt K_ETOHt ETOH ETOHe
% latex(V_ETOHt == K_ETOHt * (ETOH - ETOHe) )
% syms V_GLYCt K_GLYCt GLYC GLYCe
% latex(V_GLYCt == K_GLYCt * (GLYC - GLYCe) )
% 
% %% cofactors AXP
% 
% % AXP
% syms V_mito V_m ADP K_m_ADP PI K_m_PI % mito ATP
% latex(V_mito == V_m*ADP/(ADP+K_m_ADP)*(PI/(PI+K_m_PI)))
% % v_mito=p.mitoVmax*ADP/(ADP+p.mitoADPKm)*(PI/(PI+p.mitoPiKm)); % mito(ATP)
% 
% syms V_ATPase ATP ADP K % ATPase
% latex(V_ATPase == K * ATP/ADP)
% % v_ATPase=ATP/ADP*p.ATPaseK; % ATPase
% 
% syms V_ADK K ADP AMP ATP K_eq
% latex(V_ADK == K * (ADP^2 - (AMP*ATP)/K_eq))
% % v_ADK1=p.ADK1_k*((ADP*ADP)-(AMP*ATP)/p.ADK1_Keq);
% 
% syms V_vacPi K PI PIvac
% latex(V_vacPi == K * (PIvac - PI))
% % v_vacuolePi=p.vacuolePi_k*(p.vacuolePi_steadyStatePi-PI);
% 
% %% cofactors IXP
% syms V_Amd1 AMP K_m_AMP PI K_m_PI
% latex(V_Amd1 == (V_Amd1*AMP)/(K_m_AMP*(1+PI/K_m_PI)+AMP))
% % v_Amd1=(p.Amd1_Vmax*AMP)/(p.Amd1_K50*(1+PI/p.Amd1_Kpi)+AMP);
% syms V_Ade13Ade12 IMP K
% latex(V_Ade13Ade12 == IMP * K)
% % v_Ade13_v_Ade12=IMP*p.Ade13_Ade12_k;
% syms V_Isn1 IMP K
% latex(V_Isn1 == IMP * K)
% % v_Isn1=IMP*p.Isn1_k; 
% syms V_Pnp1 INO K
% latex(V_Pnp1 == INO * K)
% % v_Pnp1=INO*p.Pnp1_k;
% syms V_Hpt1 HYP K
% latex(V_Hpt1 == HYP * K)
% % v_Hpt1=HYP*p.Hpt1_k;
% 
% %% cofactors NADX
% syms V_mitoNADH V_m NADH K_m
% latex(V_mitoNADH == V_m * (NADH/(NADH + K_m)))
% % v_mitoNADH=p.mitoNADHVmax*(NADH/(NADH+p.mitoNADHKm));
% 
% 
% %% sink reactions
% syms V_sinkG6P V_sinkF6P V_sinkGAP V_sinkP3G V_sinkPEP V_sinkPYR V_sinkACE V_m K_m G6P F6P GAP P3G PEP PYR ACE
% latex(V_sinkG6P == V_m * (G6P/(G6P + K_m)))
% latex(V_sinkF6P == V_m * (F6P/(F6P + K_m)))
% latex(V_sinkGAP == V_m * (GAP/(GAP + K_m)))
% latex(V_sinkP3G == V_m * (P3G/(P3G + K_m)))
% latex(V_sinkPEP == V_m * (PEP/(PEP + K_m)))
% latex(V_sinkPYR == V_m * (PYR/(PYR + K_m)))
% latex(V_sinkACE == V_m * (ACE/(ACE + K_m)))
% 
% 
% %% polynomials
% syms V_m_sinkG6P V_m_sinkF6P V_m_sinkGAP V_m_sinkP3G V_m_sinkPEP V_m_sinkPYR V_m_sinkACE d
% latex(V_m_sinkG6P == 3.6854 * d.^3 -   1.4119 * d.^2 -  0.6312 * d    - 0.0043)
% latex(V_m_sinkF6P == 519.3740 * d.^6 - 447.7990 * d.^5 + 97.2843 * d.^4 + 8.0698 * d.^3 - 4.4005 * d.^2 + 0.6254 * d - 0.0078)
% latex(V_m_sinkGAP == 170.8447 * d.^6 - 113.2975 * d.^5 + 2.6494 * d.^4 + 10.2461 * d.^3 - 1.8002 * d.^2 + 0.1988 * d + 0.0012)
% latex(V_m_sinkP3G == - 0.2381 * d.^2 - 0.0210 * d - 0.0034)
% latex(V_m_sinkPEP == - 0.0637 * d.^2 - 0.0617 * d - 0.0008)
% latex(V_m_sinkPYR == - 8.4853e+03 * d.^6 + 9.4027e+03 * d.^5 - 3.8027e+03 * d.^4 + 700.5 * d.^3 - 60.26 * d.^2 + 0.711 * d - 0.0356)
% latex(V_m_sinkACE == 118.8562 * d.^6 - 352.3943 * d.^5 + 245.6092 * d.^4 - 75.2550 * d.^3 + 11.1153 * d.^2 - 1.0379 * d + 0.0119)
% % poly_sinkG6P    = 3.6854 * d.^3 -   1.4119 * d.^2 -  0.6312 * d    - 0.0043; % % CHANGED
% % poly_sinkF6P    = 519.3740 * d.^6 - 447.7990 * d.^5 + 97.2843 * d.^4 + 8.0698 * d.^3 - 4.4005 * d.^2 + 0.6254 * d - 0.0078; % % CHANGED
% % poly_sinkGAP    = 170.8447 * d.^6 - 113.2975 * d.^5 + 2.6494 * d.^4 + 10.2461 * d.^3 - 1.8002 * d.^2 + 0.1988 * d + 0.0012; % % CHANGED
% % poly_sinkP3G    = -0.2381 * d.^2 -0.0210 * d   -0.0034; % % CHANGED
% % poly_sinkPEP    = -   0.0637 * d.^2 -   0.0617 * d   -  0.0008; % % CHANGED
% % poly_sinkPYR    = - 8.4853e+03 * d.^6 + 9.4027e+03 * d.^5 - 3.8027e+03 * d.^4 + 700.5 * d.^3 - 60.26 * d.^2 + 0.711 * d - 0.0356; % % INITIAL FIT
% % poly_sinkACE    =     118.8562 * d.^6 - 352.3943 * d.^5 + 245.6092 * d.^4 - 75.2550 * d.^3 + 11.1153 * d.^2 - 1.0379 * d + 0.0119; % % CHANGED
% 
% 
% 
% %% ATPase
% %% mito(ATP)
% %% fixing vmito and vATPase to experimental values
% 
% % Now as a bypass to get it working, though unused x(88)
% if p.VmaxACE == 20
%     % v_mito
% %     RQrarioExp = [1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31];
% %     v_mito= -(v_sinkPYR * 3 + 1 * v_PDC)*1.7/2/RQrarioExp(experiment) * (ADP / (ADP+0.01));
%     v_mito_ratesSS = [0.1012    0.2248    0.4209    0.8298    1.0998    1.0958    0.7802    0.6350];
%     v_mito = v_mito_ratesSS(experiment) * (ADP / (ADP+0.01));
%     % v_ATPase
% %     v_ATPase_ratesSS = [0.2639    0.5139    0.9306    1.764    2.597    2.806    3.014    3.222];
%     v_ATPase_ratesSS = 45/40*[0.2083    0.3750    0.6528    1.2083    1.7639    1.9028    2.0417    2.1806];
%     v_ATPase = v_ATPase_ratesSS(experiment) * (ATP / (ATP+0.90));
% end
% 
% if setup.clamp10.TRE == 2
% % % % %     v_vacuolePi = zeros(size(v_vacuolePi));
%     if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))
%         v_PGM1 = zeros(size(v_PGM1));
%         v_TPS1 = zeros(size(v_TPS1));
%         v_TPS2 = zeros(size(v_TPS2));
%         v_NTH1 = zeros(size(v_NTH1));
%         v_vacuolePi = zeros(size(v_vacuolePi));
% %         v_ATPase = 0.6528 * ones(size(v_ATPase));
% %         v_mito = 0.4209 * ones(size(v_mito));
% % % % %         v_ATPase = 1.165 * ones(size(v_ATPase));
% % % % %         v_mito = 1.1026 * ones(size(v_mito));
%     end
% end
% 
% % latest changes on the ATPase and mito reactions (to make them
% % physiology-dependent for both SS and GP conditions.
% % In addition, vAmd1 re-studied here
% if setup.adjust_mito_atpase == 1
%     if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))
%         
%         % change v_mito
% %         disp('stop here');
% %         v_mito = (-v_sinkPYR * 2 + v_PDC * 1)/1.1;
%         v_mito = (-v_sinkPYR * 3 + v_PDC * 1)/1.1*0.95;
%         
%         
% % % % %             v_mito = (-v_sinkPYR * 3 + v_PDC * 1)/1.04*0.95;
%         
% %     (Vgs{j}(:,39) * 2 + 1 * Vgs{j}(:,13))/1.1
% %         (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2/[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]'
%         
% % % %         % 2021 - 01 - 11
% % % %         % NO BIOLOGICAL SENSE AT ALL, JUST TRYING P.MITOPIKM SINCE IT'S CLOSE TO 1 AND WE CAN GET IT FROM THE P ARRAY. 
% % % %         v_mito = (-v_sinkPYR * 3 + v_PDC * 1)/1.1*0.95 * p.mitoPiKm;
% % % %         v_mito = (-v_sinkPYR * 3 + v_PDC * 1)/1.1*0.95 * ADP/(ADP+p.mitoADPKm);
% 
%         % changge ATPase
%         % reference:     [0.2639    0.5139    0.9306    1.7639    2.5972    2.8056    3.0139    3.2222]
% %         v_ATPase = 0.9306;
% %         v_ATPase = 0.9306 * ATP / (ATP + p.ATPase_Katp);
%         v_ATPase = 0.9306 * ATP / ADP * p.ATPase_ratio;
% % % % %             v_ATPase = 0.9306 * ATP / ADP * 0.2707;
% 
%         % (2021/05/19) setup.final_Vmito_implementation = 1;
%         if isfield(setup,'final_Vmito_implementation')
%             if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1)&&(setup.final_Vmito_implementation == 1))
%                 v_mito=p.mitoVmax*ADP/(ADP+p.mitoADPKm)*(PI/(PI+p.mitoPiKm));
% %                 v_ATPase=ATP/ADP*p.ATPaseK;
%             end
%         end
% 
%         % (2021/03/09) unclamping IXP after pulse
%         if isfield(setup,'Katpase_fermIncrease')
%             if setup.Katpase_fermIncrease2 == 1
%                 v_ATPase = 0.9306 * ATP / ADP * p.ATPase_ratio * 3.25 * p.ATPase_ratio2;
%             
%                 if isfield(setup,'final_Vmito_implementation')
%                     if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1)&&(setup.final_Vmito_implementation == 1))
%                         v_ATPase = ATP/ADP*p.ATPaseK * 3.25 * p.ATPase_ratio2;
%                     end
%                 end
%             
%             end
%         end
%          
%         
% % % %         % 2021 - 01 - 11
% % % %         % NO BIOLOGICAL SENSE AT ALL, JUST TRYING P.MITOPIKM SINCE IT'S CLOSE TO 1 AND WE CAN GET IT FROM THE P ARRAY. 
% % % %         v_ATPase = 0.9306 * ATP / ADP * p.ATPase_ratio * p.mitoPiKm;
% 
% 
%         % % % %         v_ATPase = 0.9306 * p.mVal * ATP / ADP * p.ATPase_ratio;
% %         v_ATPase = 0.9306 * (ATP / ADP).^2 * p.ATPase_ratio;
% %         v_ATPase = 0.9306 * ATP / ADP * p.ATPase_ratio/4;
% 
%         % Change IXP
%         v_Amd1=(p.Amd1_Vmax*AMP)/(p.Amd1_K50*(1+PI/p.Amd1_Kpi)/(ATP/p.Amd1_Katp + 1)+AMP);
% %         if setup.clamp10.IXP == 1
% %             v_Amd1 = zeros(size(v_Amd1));
% %         end
% 
%         if isfield(setup,'developingIXP_20210111')
%             if setup.developingIXP_20210111 == 1
%                 % 2021 - 01 - 11 Bringing ATPase and mito closer to Y3M0 definition
%                 % during GP
%                 v_ATPase=ATP/ADP*p.ATPaseK; % ATPase
%                 v_mito=p.mitoVmax*ADP/(ADP+p.mitoADPKm)*(PI/(PI+p.mitoPiKm)); % mito(ATP)
% %                 % Y3M0
% %                 % logs
% %                 x112 = 0.4423;   
% %                 x113 = -0.2168;   
% %                 x114 = -0.0952;   
% %                 x115 = 0.0748;  
% %                 x116 = -0.1272;
% %                 x117 = 0.2991;  
% %                 x118 = 0.1285;
% %                 % pars
% %                 pAmd1_Vmax=4*10.^x112; 
% %                 pAmd1_K50=0.3*10.^x113;
% %                 pAmd1_Kpi=0.35*10.^x114;
% %                 pAde13_Ade12_k=0.05*10.^x115;
% %                 pIsn1_k=.1*10.^x116;
% %                 pPnp1_k=0.03*10.^x117;
% %                 pHpt1_k=0.02*10.^x118;
% %                 % rates
% %                v_Amd1=(pAmd1_Vmax*AMP)/(pAmd1_K50*(1+PI/pAmd1_Kpi)+AMP);
% %                v_Ade13_v_Ade12=IMP*pAde13_Ade12_k;
% %                v_Isn1=IMP*pIsn1_k; 
% %                v_Pnp1=INO*pPnp1_k;
% %                v_Hpt1=HYP*pHpt1_k;
% 
%             end
%         end
% 
%     else
%         % change v_mito
%         RQ = [1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31];
%         v_mito = (-v_sinkPYR * 3 + v_PDC * 1)/RQ(experiment)*0.95;
%         
%         % change v_ATPase
%         GAMNGAM = [0.2639    0.5139    0.9306    1.7639    2.5972    2.8056    3.0139    3.2222];
%         ATP_RATIO = [2.8280    3.2825    3.2965    3.7320    3.4540    3.3985    3.0315    2.9025];
%         ADP_RATIO = [0.7765    0.8835    0.8925    1.0175    0.8315    0.8115    0.4825    0.4795];
% %         ATPADP_RATIO = ATP_RATIO / ADP_RATIO;
% %         v_ATPase = GAMNGAM(experiment) * ATP / ADP * ATPADP_RATIO(experiment);
%         ADPATP_RATIO = ADP_RATIO / ATP_RATIO;
%         v_ATPase = GAMNGAM(experiment) * ATP / ADP * ADPATP_RATIO(experiment);
% %         v_ATPase = GAMNGAM(experiment);
%     end
%     
% elseif setup.adjust_mito_atpase == 0
% end
% 
% % check by fixing Aldolase
% if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))
%     if setup.clamp_vALDgp == 1 % if ald 1, then interpolation at t during the use of ODE
% %         disp('stop here and develop');
%         v_ALD = interp1(data.time_fluxes, data.fluxes(:,5), t, 'pchip');
%     elseif setup.clamp_vALDgp == 2 % if ald 2, final T interpolation (stop there
% %         disp('stop here and develop');
%         v_ALD = interp1(data.time_fluxes, data.fluxes(:,5), T, 'pchip');
%     end
% end
% 
% if setup.NADHrecycle_syncPYR == 1
% %     v_mitoNADH = -v_sinkPYR*3;
%     v_mitoNADH = -v_sinkPYR*3;
% end
% 
% 
% %% latest AXP-IXP implementation
% % % AXP
% % v_ADK1=p.ADK1_k*((ADP*ADP)-(AMP*ATP)/p.ADK1_Keq);
% % v_ATPase=ATP/ADP*p.ATPaseK;
% % v_mito=p.mitoVmax*ADP/(ADP+p.mitoADPKm)*(PI/(PI+p.mitoPiKm));
% % 
% % % IXP
% % v_Amd1=(p.Amd1_Vmax*AMP)/(p.Amd1_K50*(1+PI/p.Amd1_Kpi)+AMP);
% % v_Ade13_v_Ade12=IMP*p.Ade13_Ade12_k;
% % v_Isn1=IMP*p.Isn1_k; 
% % v_Pnp1=INO*p.Pnp1_k;
% % v_Hpt1=HYP*p.Hpt1_k;
% 
% % Values in FF
% 
% 
% %% temporary location
% % last edit: making Pi increase happen inside the vacoule and not in the
% % cytosol
% 
% % if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))
% %     if exist('t','var')
% %         if((t>=0)&&(t<=10))
% %             p.vacuolePi_steadyStatePi = 25;
% %         else
% %             p.vacuolePi_steadyStatePi = 10;
% %         end
% %         v_vacuolePi=p.vacuolePi_k*(p.vacuolePi_steadyStatePi-PI);
% %     elseif exist('T','var')
% %         tempP = p.vacuolePi_steadyStatePi;
% %         p.vacuolePi_steadyStatePi = 10*ones(size(PI));
% %         p.vacuolePi_steadyStatePi(1:11) = 25;
% %         v_vacuolePi=p.vacuolePi_k*(p.vacuolePi_steadyStatePi-PI);
% %     end
% % end
% 
% % % carray = [0 4 8 12 16 20 24 25 25 25 25];
% % % if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))
% % %     if exist('t','var')
% % %         if((t>=0)&&(t<=10))
% % %             p.vacuolePi_steadyStatePi = interp1([1:11],carray,t,'pchip','extrap');
% % % %             p.vacuolePi_steadyStatePi = 25;
% % %         else
% % %             p.vacuolePi_steadyStatePi = 10;
% % %         end
% % %         v_vacuolePi=p.vacuolePi_k*(p.vacuolePi_steadyStatePi-PI);
% % %     elseif exist('T','var')
% % %         tempP = p.vacuolePi_steadyStatePi;
% % %         p.vacuolePi_steadyStatePi = 10*ones(size(PI));
% % %         p.vacuolePi_steadyStatePi(1:11) = carray;
% % %         v_vacuolePi=p.vacuolePi_k*(p.vacuolePi_steadyStatePi-PI);
% % %     end
% % % end
