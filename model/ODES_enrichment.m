bmf = 2e-3; % Lin/gdw
v(1)= + v_PDC - v_ADH + v_sinkACE;% Using sinkACE (- v_ACE)
v(2)= + v_GAPDH - v_PGK;%BPG
v(3)= +v_PFK - v_ALD;%F16BP
v(4)= -v_PFK + v_PGI + v_sinkF6P + v_FRK;%F6P
v(5)= +v_GLK - v_PGI + v_sinkG6P + v_PGM1 - v_TPS1;%G6P
v(6)= -v_GLK + v_GLT + 2.*v_NTH1 + 2.*v_ATH1vac + v_glycSynthDeg + v_glycDeg; % GLCi
v(7)= + v_G3PDH - v_GAPDH + v_ADH + v_mitoNADH;%NAD
v(8)= - v_G3PDH + v_GAPDH - v_ADH - v_mitoNADH;%NADH
v(9)=  +v_ADK1 -v_GLK - v_ATPase - v_PFK +v_PGK + v_PYK - v_TPS1 + v_mito + dAXPdt.*(ATP./(ATP+ADP+AMP)) + dAXPdD.*(ATP./(ATP+ADP+AMP)) - v_FRK; %ATP
v(10)= +v_PGM - v_ENO; %P2G
v(11)= +v_PGK - v_PGM + v_sinkP3G; %P3G
v(12)= +v_ENO - v_PYK + v_sinkPEP; % %PEP
v(13)= +v_PYK - v_PDC + v_sinkPYR; % Using sinkPYR
v(14)= + v_ALD - v_GAPDH + v_TPI1 + v_sinkGAP; %GLYCERAL3P
v(15)=  -2.*v_ADK1 +v_GLK + v_ATPase + v_PFK - v_PGK - v_PYK + v_TPS1 - v_mito + dAXPdt.*(ADP./(ATP+ADP+AMP)) + dAXPdD.*(ADP./(ATP+ADP+AMP)) + v_FRK;% ADP
v(16)=  +v_ADK1 - v_Amd1 + v_Ade13_v_Ade12 + dAXPdt.*(AMP./(ATP+ADP+AMP)) + dAXPdD.*(AMP./(ATP+ADP+AMP));% AMP
v(17)= + v_ALD - v_TPI1 - v_G3PDH; % DHAP
v(18)=  + v_G3PDH - v_HOR2 - v_RHR2; % GLYC3P
v(19)= + v_HOR2 + v_RHR2 - v_GLYCEROLtransport; %GLYCEROL
v(20)= + v_ADH - v_ETOHtransport; %ETOH
v(21)= - v_PGM1 - v_UGP; %G1P
v(22) = 0; %UTP
v(23) = 0; %UDP
v(24) = + v_UGP - v_TPS1 - v_glycSynthDeg - v_glycSynth; %UDP_GLC
v(25) = + v_TPS2 - v_NTH1 - v_AGT1 - v_vacuoleT; %TRE 
v(26) = + v_TPS1 - v_TPS2; %T6P 
v(27) = - v_GAPDH + v_ATPase + v_HOR2 + v_RHR2 +2.* v_TPS1 +1.* v_TPS2 - v_mito + v_Hpt1- v_Isn1 + v_vacuolePi  + 1/2 * v_glycSynth; %PI
v(28) = +v_Amd1-v_Ade13_v_Ade12+v_Hpt1-v_Isn1; %IMP
v(29) = v_Isn1-v_Pnp1; %INO
v(30) = +v_Pnp1-v_Hpt1; %HYP
v(31) = v_ETOHtransport * Cx *bmf - ETOHec*Fout/Vbroth; % Product Ethanol
v(32) = v_GLYCEROLtransport * Cx *bmf - GLYCec*Fout/Vbroth; % Byproduct Glycerol
v(33) = + v_FRT - v_FRK; %Frci
v(34) = Fin*FRCin/Vbroth - Fout*FRCec/Vbroth - v_FRT*bmf*Cx + v_SUC; %FRCec;
v(35) = Fin*SUCin/Vbroth - Fout*SUCec/Vbroth - v_SUC; %SUCec;
v(36) = Fin*GLCin/Vbroth - Fout*GLCec/Vbroth - v_GLT*bmf*Cx + v_SUC + 2.*v_ATH1ec*bmf*Cx; %GLCec
v(37) = - Fout*TREec/Vbroth + v_AGT1*bmf*Cx - v_ATH1ec*bmf*Cx; % TREec 
v(38) = - v_ATH1vac + v_vacuoleT; % TREvac 
v(39) = (Fin - Fout); % Vbroth
v(40) = v_glycSynth - v_glycDeg; % Vbroth


%% Enrichment labeling
% labelled ratios
ACE_fL = ACE_L ./ ACE;
BPG_fL = BPG_L ./ BPG;
FBP_fL = FBP_L ./ F16BP;
F6P_fL = F6P_L ./ F6P;
G6P_fL = G6P_L ./ G6P;
GLCi_fL = GLCi_L ./ GLCi;
P2G_fL = P2G_L ./ P2G;
P3G_fL = P3G_L ./ P3G;
PEP_fL = PEP_L ./ PEP;
PYR_fL = PYR_L ./ PYR;
Glyceral3P_fL = Glyceral3P_L ./ GLYCERAL3P;
DHAP_fL = DHAP_L ./ DHAP;
Glyc3P_fL = Glyc3P_L ./ GLYC3P;
Glycerol_fL = Glycerol_L ./ GLYCEROL;
ETOH_fL = ETOH_L ./ ETOH;
G1P_fL = G1P_L ./ G1P;
UDP_GLC_fL = UDP_GLC_L ./ UDP_GLC;
TREic_fL = TREic_L ./ TRE;
T6P_fL = T6P_L ./ T6P;
ETOHec_fL = ETOHec_L ./ ETOHec;
GLYCec_fL = GLYCec_L ./ GLYCec;
FRCic_fL = FRCic_L ./ FRCi;
FRCec_fL = FRCec_L ./ FRCec;
SUCec_fL = SUCec_L ./ SUCec;
GLCec_fL = GLCec_L ./ GLCec;
TREec_fL = TREec_L ./ TREec;
TREvac_fL = TREvac_L ./ TREvac;
GLYCOGENic_fL = GLYCOGENic_L ./ Glycogen_cyt;

% 
if isnan(ACE_fL) || ACE <= 0 || ACE_L <= 0
    ACE_fL = 0;  
    ACE_L = 0;
end
if isnan(BPG_fL) || BPG<= 0 || BPG_L <= 0
    BPG_fL = 0;    
end
if isnan(FBP_fL) || F16BP<= 0 || FBP_L <= 0
    FBP_fL= 0;    
end
if isnan(F6P_fL) || F6P<= 0 || F6P_L <= 0
    F6P_fL= 0;    
end
if isnan(G6P_fL) || G6P<= 0 || G6P_L <= 0
    G6P_fL= 0;    
end
if isnan(GLCi_fL) || GLCi<= 0 || GLCi_L <= 0
    GLCi_fL= 0;    
end
if isnan(P2G_fL) || P2G<= 0 || P2G_L <= 0
    P2G_fL= 0;    
end
if isnan(P3G_fL) || P3G<= 0 || P3G_L <= 0
    P3G_fL= 0;    
end
if isnan(PEP_fL) || PEP<= 0 || PEP_L <= 0
    PEP_fL= 0;    
end
if isnan(PYR_fL) || PYR<= 0 || PYR_L <= 0
    PYR_fL= 0;    
end
if isnan(Glyceral3P_fL) || GLYCERAL3P <= 0 || Glyceral3P_L <= 0
    Glyceral3P_fL= 0;    
end
if isnan(DHAP_fL) || DHAP<= 0 || DHAP_L <= 0
    DHAP_fL= 0;    
end
if isnan(Glyc3P_fL) || GLYC3P <= 0 || Glyc3P_L <= 0
    Glyc3P_fL= 0;    
end
if isnan(Glycerol_fL) || GLYCEROL <= 0 || Glycerol_L <= 0
    Glycerol_fL= 0;    
end
if isnan(ETOH_fL) ||ETOH <= 0 || ETOH_L <= 0
    ETOH_fL= 0;    
    ETOH_L= 0; 
end
if isnan(G1P_fL) || G1P<= 0 || G1P_L <= 0
    G1P_fL= 0;    
end
if isnan(UDP_GLC_fL) || UDP_GLC<= 0 || UDP_GLC_L <= 0
    UDP_GLC_fL= 0;    
end
if isnan(TREic_fL) || TRE <= 0 || TREic_L <= 0
    TREic_fL= 0;    
end
if isnan(T6P_fL) || T6P<= 0 || T6P_L <= 0
    T6P_fL= 0;    
end
if isnan(ETOHec_fL) || ETOHec<= 0 || ETOHec_L <= 0
    ETOHec_fL= 0;  
end
if isnan(GLYCec_fL) || GLYCec<= 0 || GLYCec_L <= 0
    GLYCec_fL= 0;    
end
if isnan(FRCic_fL) || FRCi <= 0 || FRCic_L <= 0
    FRCic_fL= 0;    
end
if isnan(FRCec_fL) || FRCec<= 0 || FRCec_L <= 0
    FRCec_fL= 0;    
end
if isnan(SUCec_fL) || SUCec<= 0 || SUCec_L <= 0
    SUCec_fL= 0;    
end
if isnan(GLCec_fL) || GLCec<= 0 || GLCec_L <= 0
    GLCec_fL= 0;    
end
if isnan(TREec_fL) || TREec<= 0 || TREec_L <= 0
    TREec_fL= 0;    
end
if isnan(TREvac_fL) || TREvac<= 0 || TREvac_L <= 0
    TREvac_fL= 0;
end


% new fluxes: fwd and bwd for ADH, ETOH_transport, AGT1, FRT, FRK
v(41) = + PYR_fL .* v_PDC - ACE_fL .* v_ADH + ACE_fL .* v_sinkACE ; % ACE#1 :  % missing sink  + ETOH#1 * r_41_v_ADH_bwd
v(42) = + Glyceral3P_fL .* v_GAPDH - BPG_fL .* v_PGK; % BPG#1 :  
v(43) = + F6P_fL .* v_PFK - FBP_fL .* v_ALD ; % FBP#1 :  
v(44) = + G6P_fL .* v_PGI - F6P_fL .* v_PFK + FRCic_fL .* v_FRK + F6P_fL .* v_sinkF6P; % F6P#1 % missing sink    - F6P#1 * r_44_v_FRK_bwd
v(45) = + GLCi_fL .* v_GLK - G6P_fL .* v_PGI + G6P_fL .* v_PGM1 - G6P_fL .* v_TPS1 + G6P_fL .* v_sinkG6P; % G6P#1 :  % missing sink
v(46) = + GLCec_fL .* v_GLT - GLCi_fL .* v_GLK + TREic_fL .* 2.* v_NTH1 + TREvac_fL .* 2.* v_ATH1vac  + v_glycSynthDeg * GLYCOGENic_fL + v_glycDeg * GLYCOGENic_fL; % GLCic#1 :  
v(47) = + P3G_fL .* v_PGM - P2G_fL .* v_ENO; % P2G#1 :
v(48) = + BPG_fL .* v_PGK - P3G_fL .* v_PGM + P3G_fL .* v_sinkP3G; % P3G#1 :  %missing sink
v(49) = + P2G_fL .* v_ENO - PEP_fL .* v_PYK + PEP_fL .* v_sinkPEP; % PEP#1 :  %missing sink
v(50) = + PEP_fL .* v_PYK - PYR_fL .* v_PDC + PYR_fL .* v_sinkPYR; % PYR#1 : %missing sink 
v(51) = + FBP_fL .* v_ALD + DHAP_fL .* v_TPI1 - Glyceral3P_fL .* v_GAPDH + Glyceral3P_fL .* v_sinkGAP; % Glyceral3P#1 :  %missing sink
v(52) = + FBP_fL .* v_ALD - DHAP_fL .* v_G3PDH - DHAP_fL .* v_TPI1; % DHAP#1 :  
v(53) = + DHAP_fL .* v_G3PDH - Glyc3P_fL .* v_HOR2 - Glyc3P_fL .* v_RHR2 ; %Glyc3P#1 :  
v(54) = + Glyc3P_fL .* v_HOR2 + Glyc3P_fL .* v_RHR2 - Glycerol_fL .* v_GLYCEROLtransport; % Glycerol#1 :  
v(55) = - ETOH_fL .* v_ETOHtransport + ACE_fL .* v_ADH ;% ETOH#1 : - ETOH#1 * r_41_v_ADH_bwd + ETOHec#1 * r_23_v_ETOHtransport_bwd 
v(56) = - G6P_fL .* v_PGM1 - G1P_fL .* v_UGP; % G1P#1 :  
v(57) = + G1P_fL .* v_UGP - UDP_GLC_fL .* (v_TPS1 + v_glycSynthDeg + v_glycSynth); % UDP_GLC#1 :
v(58) =  + T6P_fL .* v_TPS2 - TREic_fL .* v_NTH1 - TREic_fL .* v_AGT1 - TREic_fL .* v_vacuoleT ; % TREic#1 : + TREec#1 * r_47_v_AGT1_bwd 
v(59) =  - T6P_fL .* v_TPS2 + (6 .* G6P_fL + 6 .* UDP_GLC_fL) ./ 12 .* v_TPS1; % T6P#1 : 
v(60) =  + ETOH_fL .* v_ETOHtransport .*Cx.*bmf - ETOHec_fL .* ETOHec.*Fout./Vbroth ; % ETOH#1 extracellular  - ETOHec#1 * r_23_v_ETOHtransport_bwd
v(61) = Glycerol_fL .* v_GLYCEROLtransport .*Cx.*bmf - GLYCec_fL .* GLYCec.*Fout./Vbroth ;  % GLYCec#1 extracellular 
v(62) =  + FRCec_fL .* v_FRT - FRCic_fL .* v_FRK; % FRUCic#1 : + F6P#1 * r_44_v_FRK_bwd - FRUCic#1 * r_43_v_FRT_bwd
v(63) =  - FRCec_fL .* v_FRT .*Cx.*bmf + SUCec_fL .* v_SUC + FRCin_fL .* Fin.*FRCin./Vbroth - FRCec_fL .* Fout.*FRCec./Vbroth; % FRUCec#1 : + FRUCic#1 * r_43_v_FRT_bwd
v(64) =  SUCin_fL .* Fin.*SUCin./Vbroth - SUCec_fL .* Fout.*SUCec./Vbroth - SUCec_fL .* v_SUC; %SUCec
v(65) =  - GLCec_fL .* v_GLT .*Cx.*bmf + TREec_fL .* 2.* v_ATH1ec .*Cx.*bmf + SUCec_fL .* v_SUC + GLCin_fL .* Fin.*GLCin./Vbroth - GLCec_fL .* Fout.*GLCec./Vbroth; % GLCec#1 :  % (+ Feed_fL .* vFeed;)
v(66) =  - TREec_fL .* v_ATH1ec .*Cx.*bmf + TREic_fL .* v_AGT1 .*Cx.*bmf + TREin_fL .* Fin.*TREin./Vbroth - TREec_fL .* Fout.*TREec./Vbroth ; % TREec#1 : - TREec#1 * r_47_v_AGT1_bwd
v(67) =  - TREvac_fL .* v_ATH1vac + TREic_fL .* v_vacuoleT; % TREvac#1 : 
v(68) = + UDP_GLC_fL .* v_glycSynth - GLYCOGENic_fL .* v_glycDeg;


%% Optional: Setting metabolite dcdts to zero - select by comment/uncommenting
NI = [];

% TREic
% NI = 25;

% % Only upper glycolysis until FBP
% NI = [1:2, 7:32];

% % Upper glyc + GAP&DHAP
% NI = [1:2, 7:13, 15:16, 18:32]; 

% % Upper glyc + GAP&DHAP + BPG
% NI = [1, 7:13, 15:16, 18:32];

% % Upper glyc + lower until PEP
% NI = [1, 7:9, 13, 15:16, 18:32];

% % Upper glyc + lower until PYR
% NI = [1, 7:9, 15:16, 18:32];

% % Upper glyc + lower until PYR, including acetate
% NI = [7:9, 15:16, 18:32];

%% Including glycerol and ethanol
% % Only upper glycolysis until FBP
% NI = [1:2, 7:18, 21:30];

% % Upper glyc + GAP&DHAP
% NI = [1:2, 7:13, 15:16, 18, 21:30]; 

% % Upper glyc + GAP&DHAP + BPG
% NI = [1, 7:13, 15:16, 18, 21:30];

% % Upper glyc + lower until PEP
% NI = [1, 7:9, 13, 15:16, 18, 21:30];

% % Upper glyc + lower until PYR
% NI = [1, 7:9, 15:16, 18, 21:30];

% % Upper glyc + lower until PYR, including acetate
% NI = [7:9, 15:16, 18, 21:30];

%% Including glycerol and ethanol and glyc3p
% % Only upper glycolysis until FBP
% NI = [1:2, 7:17, 21:30];

% % Upper glyc + GAP&DHAP
% NI = [1:2, 7:13, 15:16, 21:30]; 

% % Upper glyc + GAP&DHAP + BPG
% NI = [1, 7:13, 15:16, 21:30];

% % % Upper glyc + lower until PEP
% NI = [1, 7:9, 13, 15:16, 21:30];

% Upper glyc + lower until PYR
% NI = [1, 7:9, 15:16, 21:30];

% % Upper glyc + lower until PYR, including acetate
%NI = [7:9, 15:16, 21:30];

%% Including glycerol and ethanol and glyc3p and acetate
% % Only upper glycolysis until FBP
%NI = [2, 7:17, 21:30];

% % % Upper glyc + GAP&DHAP
% NI = [2, 7:13, 15:16, 21:30]; 

% % Upper glyc + GAP&DHAP + BPG
%NI = [7:13, 15:16, 21:30];

% % Upper glyc + lower until PEP
%NI = [7:9, 13, 15:16, 21:30];

% % Upper glyc + lower until PYR, including acetate
%NI = [7:9, 15:16, 21:30];

%% everything except pyruvate 
% %NI = 13;
% %% Everything excep glycerol, ethanol
% NI = [19, 20, 31, 32]; 
% 
% %% Everything excep glycerol, ethanol, trehalose
% NI = [19, 20, 31, 32, 25]; 
v(NI) = 0; 