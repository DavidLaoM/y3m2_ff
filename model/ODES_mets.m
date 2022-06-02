if setup.GPdataset.GP400WT        == 1
    GlucoseFeed400s;
elseif setup.GPdataset.GP1800WT 	== 1
    GlucoseFeed1800s;
elseif setup.GPdataset.GP400M     == 1
    GlucoseFeed400s;
elseif setup.GPdataset.Sucrose == 1
    SucroseFeed400s;
elseif setup.GPdataset.Fructose == 1
    FructoseFeed400s;
end

bmf = 2e-3; % Lin./gdw
v(:,1)= + v_PDC - v_ADH + v_sinkACE;% Using sinkACE (- v_ACE)
v(:,2)= + v_GAPDH - v_PGK;%BPG
v(:,3)= +v_PFK - v_ALD;%F16BP
v(:,4)= -v_PFK + v_PGI + v_sinkF6P + v_FRK;%F6P
v(:,5)= +v_GLK - v_PGI + v_sinkG6P;%G6P                                       %TRE cycle off here:
v(:,6)= -v_GLK + v_GLT + 2.*v_NTH1 + 2.*v_ATH1vac;  %GLCi
v(:,7)= + v_G3PDH - v_GAPDH + v_ADH + v_mitoNADH;%NAD
v(:,8)= - v_G3PDH + v_GAPDH - v_ADH - v_mitoNADH;%NADH
v(:,9)=  +v_ADK1 -v_GLK - v_ATPase - v_PFK +v_PGK + v_PYK - v_TPS1 + v_mito + dAXPdt.*(ATP./(ATP+ADP+AMP)) + dAXPdD.*(ATP./(ATP+ADP+AMP)); %ATP
v(:,10)= +v_PGM - v_ENO; %P2G
v(:,11)= +v_PGK - v_PGM + v_sinkP3G; %P3G
v(:,12)= +v_ENO - v_PYK + v_sinkPEP; % %PEP
v(:,13)= +v_PYK - v_PDC + v_sinkPYR; % Using sinkPYR
v(:,14)= + v_ALD - v_GAPDH + v_TPI1 + v_sinkGAP; %GLYCERAL3P
v(:,15)=  -2.*v_ADK1 +v_GLK + v_ATPase + v_PFK - v_PGK - v_PYK + v_TPS1 - v_mito + dAXPdt.*(ADP./(ATP+ADP+AMP)) + dAXPdD.*(ADP./(ATP+ADP+AMP));% ADP
v(:,16)=  +v_ADK1 - v_Amd1 + v_Ade13_v_Ade12 + dAXPdt.*(AMP./(ATP+ADP+AMP)) + dAXPdD.*(AMP./(ATP+ADP+AMP));% AMP
v(:,17)= + v_ALD - v_TPI1 - v_G3PDH; % DHAP
v(:,18)=  + v_G3PDH - v_HOR2 - v_RHR2; % GLYC3P
v(:,19)= + v_HOR2 + v_RHR2 - v_GLYCEROLtransport; %GLYCEROL
v(:,20)= + v_ADH - v_ETOHtransport; %ETOH
v(:,21)= - v_PGM1 - v_UGP; %G1P
v(:,22) = 0; %UTP
v(:,23) = 0; %UDP
v(:,24) = + v_UGP - v_TPS1; %UDP_GLC
v(:,25) = + v_TPS2 - v_NTH1 + v_AGT1 - v_vacuoleT; %TRE 
v(:,26) = + v_TPS1 - v_TPS2; %T6P 
v(:,27) = - v_GAPDH + v_ATPase + v_HOR2 + v_RHR2 +2.* v_TPS1 +1.* v_TPS2 - v_mito + v_Hpt1- v_Isn1 + v_vacuolePi; %PI
v(:,28) = +v_Amd1-v_Ade13_v_Ade12+v_Hpt1-v_Isn1; %IMP
v(:,29) = v_Isn1-v_Pnp1; %INO
v(:,30) = +v_Pnp1-v_Hpt1; %HYP
v(:,31) = v_ETOHtransport .* Cx .*bmf - ETOHec.*Fout./Vbroth; % Product Ethanol
v(:,32) = v_GLYCEROLtransport .* Cx .*bmf - GLYCec.*Fout./Vbroth; % Byproduct Glycerol

v(:,33) = + v_FRT - v_FRK; %Frci
v(:,34) = Fin.*FRCin./Vbroth - Fout.*FRCec./Vbroth - v_FRT.*bmf.*Cx + v_SUC; %FRCec;
v(:,35) = Fin.*SUCin./Vbroth - Fout.*SUCec./Vbroth - v_SUC; %SUCec;
v(:,36) = Fin.*GLCin./Vbroth - Fout.*GLCec./Vbroth - v_GLT.*bmf.*Cx + v_SUC + 2.*v_ATH1ec; %GLCec
v(:,37) = - Fout.*TREec./Vbroth - v_AGT1.*bmf.*Cx - v_ATH1ec; % TREec 
v(:,38) = - v_ATH1vac + v_vacuoleT; % TREvac
v(:,39) = (Fin - Fout);
