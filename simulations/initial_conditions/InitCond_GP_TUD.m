%% From Suarez Mendez appendix 2.2 - D = 0.101 h-1
z = 2; % 1gdw = 2ml

TREec = 0;
TREvac = 145.94; % 90% of the measured trehalose
% TREvac = 16.22; % 10% of the measured trehalose

% Extracellular Concentrations
GLCec = 0.183; % mM
FRCec = 0; % mM
SUCec = 0; % mM
ACEec = 3.58; % mM
ETOHec = 0;
GLYCec =0; 
FRCi = 0;
% Glycolysis - umol/gdw
GLCi = 5.57;
G6P = 5.42; 
F6P = 1.5; 
M6P = 1.76;
F16BP = 0.23;
DHAP = 0.54;
GAP = 0.032;
G3P = 0.207;
P3G = 4.99; 
P2G = 0.587;
PEP = 2.41;
PYR = 2.20; 

GLYCERAL3P = GAP;
GLYC3P = G3P;
% Storage Carbohydrate Pools
T6P = 0.374;
TRE = 16.22; % 10% of the measured trehalose
% TRE = 145.94; % 90% of the measured trehalose
Glycogen = 593;
G1P = 0.36;
UDP_GLC = 3.85;

% PPP
CIT = 6.96;
iCIT= 0.40;
AKG = 1.82;
SUCC = 0.96;
FUM = 0.67;
MAL = 3.01;

% Others
AcCoA = 0.25;
CoA = 0.12;

% Cofactors
NADH = 0.47;
NAD = 1.83;
NADPH = 1.88;
NADP = 0.41;
FAD = 0.13;

% Nucleotides 
AMP = 0.55;
ADP = 2.15;
ATP = 9.70;
AXP = AMP + ADP + ATP; 

UMP = 0.11;
UDP = 0.74;
UTP = 1.84; 
UXP = UMP + UDP + UTP; 

Vbroth = 3.894; %L

%% notes
% GLYCERATE13BP=BPG;
%% Unknowns

ACE=0.04;               % Not available. Left the same way.
BPG=0;                  % Not available. Left the same way.
GLYCEROL=0.1;           % Not available. Left the same way.  %glycorol [cytoplasm] (mmol/L)
ETOH=0.1;                % Not available. Left the same way. %25;

PI=20;  % could be higher (15-20)                % Not available. Left the same way.
IMP=0.1;%0.5; %IMP      % Not available. Left the same way.
INO=0.1;%0.5; %inosine  % Not available. Left the same way.
HYP=1.5; %1 %hypoxyanthine  % Not available. Left the same way.



IC=[ACE;
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
ETOHec*z;
GLYCec*z;
FRCi;
FRCec*z;
SUCec*z;
GLCec*z;
TREec*z;
TREvac;
Vbroth*z;]/z;

Glycogen_cyt = 100;
IC(40) = Glycogen_cyt;

f.GLYCEROL_e=GLYCec; %glycerol_out [extracellular] (mmol/L)
f.ETOH_e=ETOHec; %ethanol_out [extracellular] (mmol/L)

clear ACE BPG F16BP F6P G6P GLCi NAD NADH ATP P2G P3G PEP PYR GLYCERAL3P ADP
clear AMP DHAP GLYC3P GLYCEROL ETOH G1P UTP UDP UDP_GLC TRE T6P PI IMP INO HYP
clear ETOHec GLYCec FRCi SUCec FRCec GLCec TREec TREvac
clear AcCoA ACEec AKG AXP CIT CoA FAD FUM G3P GAP Glycogen M6P MAL NADP NADPH SUCC UMP UXP