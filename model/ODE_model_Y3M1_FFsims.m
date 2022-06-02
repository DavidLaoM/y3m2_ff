function v = ODE_model_Y3M1_FFsims(t,IC,p,f,d,setup, data)
% system of ODE equations in the Y3M1 model. The clamps affect the
% calculation of reaction rates in rateEquations_Y3M1.

% recall metabolite concentrations
IC = real(IC);
ACE=IC(1);
BPG=IC(2);
F16BP=IC(3);
F6P=IC(4);
G6P=IC(5);
GLCi=IC(6);
NAD=IC(7);
NADH=IC(8);
ATP=IC(9); 
P2G=IC(10);
P3G=IC(11);
PEP=IC(12);
PYR=IC(13);
GLYCERAL3P=IC(14);
ADP=IC(15);
AMP=IC(16); 
DHAP=IC(17);
GLYC3P=IC(18);
GLYCEROL=IC(19);
ETOH=IC(20);
G1P=IC(21);
UTP=IC(22);
UDP=IC(23);
UDP_GLC=IC(24);
TRE=IC(25);
T6P=IC(26);
PI=IC(27);
IMP=IC(28);
INO=IC(29);
HYP=IC(30);
ETOHec=IC(31);
GLYCec=IC(32);
FRCi = IC(33); 
FRCec = IC(34);
SUCec = IC(35);
GLCec = IC(36);
TREec = IC(37);
TREvac = IC(38);
Vbroth = IC(39);
Glycogen_cyt = IC(40);

% clamps
if setup.stage == 1
    if setup.GPdataset.GP400WT == 1
        UDP_GLC_exp = data.FF01.metabolites.ICUDPG.conc(1);
    elseif setup.GPdataset.GP1800WT == 1
        UDP_GLC_exp = data.FF03.metabolites.ICUDPG.conc(1);
    elseif setup.GPdataset.GP400M == 1
        UDP_GLC_exp = data.FF04.metabolites.ICUDPG.conc(1);
    end
elseif setup.stage == 2
    if setup.GPdataset.GP400WT == 1
        UDP_GLC_exp = interp1(data.FF01.metabolites.ICUDPG.time,data.FF01.metabolites.ICUDPG.conc,t,'pchip','extrap');
    elseif setup.GPdataset.GP1800WT == 1
        UDP_GLC_exp = interp1(data.FF03.metabolites.ICUDPG.time,data.FF03.metabolites.ICUDPG.conc,t,'pchip','extrap');
    elseif setup.GPdataset.GP400M == 1
        UDP_GLC_exp = interp1(data.FF04.metabolites.ICUDPG.time,data.FF04.metabolites.ICUDPG.conc,t,'pchip','extrap');
    end
end

% extracellular glucose and experimental setup conditions
f.GLCo = GLCec;
if setup.GPdataset.GP400WT        == 1
    GlucoseFeed400s;
elseif setup.GPdataset.GP1800WT 	== 1
    GlucoseFeed1800s;
elseif setup.GPdataset.GP400M     == 1
    GlucoseFeed400s;
end

% Clamping option for extracellular glucose concentration
if((isfield(setup,'clamp_GLCec'))&&(setup.clamp_GLCec == 1)&&(setup.stage == 2))
    temp = data.FF01.metabolites.ECglucose;
    GLCec = interp1(temp.time(1:15), temp.conc(1:15), t, 'pchip');
end

% reaction rates
rateEquations_Y3M1;

% clamping rates option
if setup.clamp.TRE == 1 % trehalose cycle
    v_PGM1 = zeros(size(v_PGM1));
    v_TPS1 = zeros(size(v_TPS1));
    v_TPS2 = zeros(size(v_TPS2));
    v_NTH1 = zeros(size(v_NTH1));
    v_UGP = zeros(size(v_UGP));
    v_ATH1ec = zeros(size(v_ATH1ec));
    v_ATH1vac = zeros(size(v_ATH1vac));
    v_AGT1 = zeros(size(v_AGT1));
    v_vacuoleT = zeros(size(v_vacuoleT));
end
if setup.clamp.IXP == 1 % ixp
    v_Amd1 = zeros(size(v_Amd1));
    v_Ade13_v_Ade12 = zeros(size(v_Ade13_v_Ade12));
    v_Isn1 = zeros(size(v_Isn1));
    v_Pnp1 = zeros(size(v_Pnp1));
    v_Hpt1 = zeros(size(v_Hpt1));
    dAXPdt = 0; % meanwhile developing AXP and IXP
    dAXPdD = 0; % meanwhile developing AXP and IXP
end 
if((isfield(setup.GPdataset,'GP1800WT'))&&(setup.GPdataset.GP1800WT == 1)) % specific of GP1800WT setup (missing experimental data)
    if isreal(v_PDC) == 0
        v_PDC = zeros(size(v_PDC));
    end
end

% system of ODEs
ODES

v = v';
if((isfield(setup.GPdataset,'GP1800WT'))&&(setup.GPdataset.GP1800WT == 1)) % specific of GP1800WT setup (missing experimental data)
    if isreal(v) == 0
        fprintf('t=%f, vreal? %f. \n',t,isreal(v))
    end
end

end