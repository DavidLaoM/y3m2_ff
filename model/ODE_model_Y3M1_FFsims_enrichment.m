    function v = ODE_model_Y3M1_FFsims_enrichment(t,IC,p,f,d,setup, data)
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

% added enrichment IC values
ACE_L = IC(41);
BPG_L = IC(42);
FBP_L = IC(43);
F6P_L = IC(44);
G6P_L = IC(45);
GLCi_L = IC(46);
P2G_L = IC(47);
P3G_L = IC(48);
PEP_L = IC(49);
PYR_L = IC(50);
Glyceral3P_L = IC(51);
DHAP_L = IC(52);
Glyc3P_L = IC(53);
Glycerol_L = IC(54);
ETOH_L = IC(55);
G1P_L = IC(56);
UDP_GLC_L = IC(57);
TREic_L = IC(58);
T6P_L = IC(59);
ETOHec_L = IC(60);
GLYCec_L = IC(61);
FRCic_L = IC(62);
FRCec_L = IC(63);
SUCec_L = IC(64);
GLCec_L = IC(65);
TREec_L = IC(66);
TREvac_L = IC(67);
GLYCOGENic_L = IC(68);

% extracellular glucose and experimental setup options
f.GLCo = GLCec;
if setup.GPdataset.GP400WT        == 1
    GlucoseFeed400s;
elseif setup.GPdataset.GP1800WT 	== 1
    GlucoseFeed1800s;
elseif setup.GPdataset.GP400M     == 1
    GlucoseFeed400s;
end

% clamps
if setup.stage == 1
    UDP_GLC_exp = data.FF01.metabolites.ICUDPG.conc(1);
elseif setup.stage == 2
    UDP_GLC_exp = interp1(data.FF01.metabolites.ICUDPG.time,data.FF01.metabolites.ICUDPG.conc,t,'pchip','extrap');
end

% clamping of labelled GLCec
if((isfield(setup,'clamp_enrichment_GLCec'))&&(setup.clamp_enrichment_GLCec == 1))
    temp_1 = setup.clamp_enrichment_GLCec_data.metabolites{65};
    temp_fraction = interp1(temp_1.time, temp_1.fraction, t,'pchip');
    GLCec_L = GLCec * temp_fraction;
end

% calculate reaction rates
rateEquations_Y3M1_enrichment;

% clamps
if setup.clamp.TRE == 1 % Trehalose cycle
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
if setup.clamp.IXP == 1 % IXP cycle
    v_Amd1 = zeros(size(v_Amd1));
    v_Ade13_v_Ade12 = zeros(size(v_Ade13_v_Ade12));
    v_Isn1 = zeros(size(v_Isn1));
    v_Pnp1 = zeros(size(v_Pnp1));
    v_Hpt1 = zeros(size(v_Hpt1));
    dAXPdt = 0; % meanwhile developing AXP and IXP
    dAXPdD = 0; % meanwhile developing AXP and IXP
end

% system of ODEs
ODES_enrichment

v = v';
end