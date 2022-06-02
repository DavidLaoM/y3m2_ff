%% IC values recall
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
IMP=IC(8);
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


%% Labelled mets 
ACE_L = 0.01 * ACE;
BPG_L = 0.01 * BPG;
FBP_L = 0.01 * F16BP;
F6P_L = 0.01 * F6P;
G6P_L = 0.01 * G6P;
GLCi_L = 0.01 * GLCi;
P2G_L = 0.01 * P2G;
P3G_L = 0.01 * P3G;
PEP_L = 0.01 * PEP;
PYR_L = 0.01 * PYR;
Glyceral3P_L = 0.01 * GLYCERAL3P;
DHAP_L = 0.01 *DHAP ;
Glyc3P_L = 0.01 * GLYC3P ;
Glycerol_L = 0.01 *GLYCEROL ;
ETOH_L = 0.01 * ETOH;
G1P_L = 0.01 *G1P ;
UDP_GLC_L = 0.01 *UDP_GLC ;
TREic_L = 0.01 *TRE ;
T6P_L = 0.01 *T6P ;
ETOHec_L = 0.01 * ETOHec;
GLYCec_L = 0.01 * GLYCec;
FRCic_L = 0.01 * FRCi;
FRCec_L = 0.01 *FRCec ;
SUCec_L = 0.01 *SUCec ;
GLCec_L = 0.01 * GLCec;
TREec_L = 0.01 *TREec ;
TREvac_L = 0.01 *TREvac ;


%% merging them in IC-vector
clear IC
IC=[ACE; % first non enriched data
BPG;
F16BP;
F6P;
G6P;
GLCi;
NAD;
NADH;
ATP;
P2G;
P3G;
PEP;
PYR;
GLYCERAL3P;
ADP;
AMP;
DHAP;
GLYC3P;
GLYCEROL;
ETOH;
G1P;
UTP;
UDP;
UDP_GLC;
TRE;
T6P;
PI;
IMP;
INO;
HYP;
ETOHec;
GLYCec;
FRCi;
FRCec;
SUCec;
GLCec;
TREec;
TREvac;
Vbroth;
Glycogen_cyt;
ACE_L; % enriched from here onwards
BPG_L;
FBP_L;
F6P_L;
G6P_L;
GLCi_L;
P2G_L;
P3G_L;
PEP_L;
PYR_L;
Glyceral3P_L;
DHAP_L;
Glyc3P_L;
Glycerol_L;
ETOH_L;
G1P_L;
UDP_GLC_L;
TREic_L;
T6P_L;
ETOHec_L;
GLYCec_L;
FRCic_L;
FRCec_L;
SUCec_L;
GLCec_L;
TREec_L;
TREvac_L;
0];


f.GLYCEROL_e=GLYCec; %glycerol_out [extracellular] (mmol/L)
f.ETOH_e=ETOHec; %ethanol_out [extracellular] (mmol/L)

clear ACE BPG F16BP F6P G6P GLCi NAD NADH ATP P2G P3G PEP PYR GLYCERAL3P ADP
clear AMP DHAP GLYC3P GLYCEROL ETOH G1P UTP UDP UDP_GLC TRE T6P PI IMP INO HYP
clear ETOHec GLYCec FRCi SUCec FRCec GLCec TREec TREvac
clear AcCoA ACEec AKG AXP CIT CoA FAD FUM G3P GAP Glycogen M6P MAL NADP NADPH SUCC UMP UXP

clear ACE_L BPG_L FBP_L F6P_L G6P_L GLCi_L P2G_L P3G_L PEP_L PYR_L Glyceral3P_L DHAP_L
clear Glyc3P_L Glycerol_L ETOH_L G1P_L UDP_GLC_L TREic_L T6P_L ETOHec_L
clear GLYCec_L FRCic_L FRCec_L SUCec_L GLCec_L TREec_L TREvac_L