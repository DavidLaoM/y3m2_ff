ACE=0.04;               % Not available. Left the same way.
BPG=0;                  % Not available. Left the same way.
%convert to mmol/L assuming 1 gdw = 2 mL (divide excel IC concentration per
%2).
F16BP   = 0.015170432039596; %0.18;
F6P     = 0.173387606382047; %0.46;
G6P     = 0.61023414420927; %1.7;
GLCi    = 0.927858622293273/10; % Assuming 10% of total measured GLCi %= 0.2;               % Not available. Left the same way. 
NAD=1.49;               % Not available. Left the same way.
NADH=0.01;              % Not available. Left the same way.
ATP     = 5.12685120427421; %4.43.;
P2G     = 0.427631746363227; %0.1;
P3G     = 3.72351228834017; %3.1;
PEP     = 2.03246228675531; %1.4;
PYR     = 1.56384844404728; %1;
GLYCERAL3P = 0.0141849231126479; %0.007; %glyceraldehyde 3-phosphate [cytoplasm] (mmol/L)
ADP     = 5.27515835489072; %1.06.;
AMP     = 2.26632223119646; %0.0693.;
DHAP    = 0.0435910821743525; %1; %glycerone phosphate [cytoplasm] (mmol/L)
GLYC3P  = 0.00040193212621191; %0.2; %sn-glycerol 3-phosphate [cytoplasm] (mmol/L)
GLYCEROL=0.1;           % Not available. Left the same way.  %glycorol [cytoplasm] (mmol/L)
ETOH=10;                % Not available. Left the same way. %25;
G1P     = 0.0772415779162459; %0.11; %D-glucose 1-phosphate [cytoplasm] (mmol/L)
UTP=0.649;              % Not available. Left the same way. %UTP [cytoplasm] (mmol/L)
UDP=0.281;              % Not available. Left the same way. %UDP [cytoplasm] (mmol/L)
UDP_GLC = 1.18782065341151; %0.07; %UPD-D-glucose [cytoplasm] (mmol/L)
TRE     = 240.446632812843; %56; %alpha,alpha-trehalose [cytoplasm] (mmol/L)
T6P     = 0.356775196005687; %0.11; %alpha,aplha-trehalose 6-phosphate [cytoplasm] (mmol/L)
PI=20;                  % Not available. Left the same way.
IMP=0.1;%0.5; %IMP      % Not available. Left the same way.
INO=0.1;%0.5; %inosine  % Not available. Left the same way.
HYP=1.5; %1 %hypoxyanthine  % Not available. Left the same way.
FRCi = 0;
FRCec = 0;
SUCec = 0;
GLCec = 0.169;
Vbroth = 3.9;
ETOHec = 0;
GLYCec = 0;

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
G1P
UTP;
UDP;
UDP_GLC;
TRE;
T6P;
PI;
IMP;
INO;
HYP;
FRCi;
FRCec;
SUCec;
GLCec; 
Vbroth];IC=[ACE;
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
Vbroth*z]/z;

%fixed states           % Not available. Left the same way.
f.CO2=1;
% % % % % f.GLCo=0.0943666794; %0.1;%110;
f.GLY=0.15;
f.Glyc=0;
f.SUCC=0;
f.Trh=0;
f.X=0;
f.GLYCEROL_e=0; %glycerol_out [extracellular] (mmol/L)
f.ETOH_e=0; %ethanol_out [extracellular] (mmol/L)
f.TRE_e=0; %alpha,alpha-trehalose_out [extracellular] (mmol/L)

clear ACE BPG F16BP F6P G6P GLCi NAD NADH ATP P2G P3G PEP PYR GLYCERAL3P ADP
clear AMP DHAP GLYC3P GLYCEROL ETOH G1P UTP UDP UDP_GLC TRE T6P PI IMP INO HYP
clear ETOHec GLYCec FRCi