z = 2;
%convert to mmol/L assuming 1 gdw = 2 mL (divide excel IC concentration per
%2).
F16BP   = 0.0361236129771379; %0.18;
F6P     = 0.0630922905446206; %0.46;
G6P     = 0.233205459132109; %1.7;
GLCi    = 0.290871750979068/10; % Assuming 10% of total measured GLCi %= 0.2;               % Not available. Left the same way. 

ATP     = 3.3516417428141; %4.43.;
P2G     = 0.342165476917415; %0.1;
P3G     = 2.83185839233627; %3.1;
PEP     = 1.61061750787239; %1.4;
PYR     = 0.95226705421595; %1;
GLYCERAL3P = 0.0114915085230956; %0.007; %glyceraldehyde 3-phosphate [cytoplasm] (mmol/L)
ADP     = 2.77627477977165; %1.06.;
AMP     = 1.21639639826066; %0.0693.;
DHAP    = 0.0230940866590949; %1; %glycerone phosphate [cytoplasm] (mmol/L)
GLYC3P  = 0.00718396836231064; %0.2; %sn-glycerol 3-phosphate [cytoplasm] (mmol/L)

G1P     = 0.0970148151499796; %0.11; %D-glucose 1-phosphate [cytoplasm] (mmol/L)

UDP_GLC = 0.818534179166591; %0.07; %UPD-D-glucose [cytoplasm] (mmol/L)
TRE     = 44.937152473174; %56; %alpha,alpha-trehalose [cytoplasm] (mmol/L)
T6P     = 0.113893674303651; %0.11; %alpha,aplha-trehalose 6-phosphate [cytoplasm] (mmol/L)

FRCi = 0;

GLCec = 0.493;
Vbroth = 3.93; 

%% Unknowns
NAD=1.49;               % Not available. Left the same way.
NADH=0.01;              % Not available. Left the same way.
GLYCEROL=0.1;           % Not available. Left the same way.  %glycorol [cytoplasm] (mmol/L)
ETOH=10;                % Not available. Left the same way. %25;
UTP=0.649;              % Not available. Left the same way. %UTP [cytoplasm] (mmol/L)
UDP=0.281;              % Not available. Left the same way. %UDP [cytoplasm] (mmol/L)
PI=20;                  % Not available. Left the same way.
IMP=0.1;%0.5; %IMP      % Not available. Left the same way.
INO=0.1;%0.5; %inosine  % Not available. Left the same way.
HYP=1.5; %1 %hypoxyanthine  % Not available. Left the same way.
ACE=0.04;               % Not available. Left the same way.
BPG=0;                  % Not available. Left the same way.

ETOHec = 0;
GLYCec = 0;
FRCec = 0;
SUCec = 0;


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
Vbroth*z]/z;


f.GLYCEROL_e=GLYCec; %glycerol_out [extracellular] (mmol/L)
f.ETOH_e=ETOHec; %ethanol_out [extracellular] (mmol/L)

clear ACE BPG F16BP F6P G6P GLCi NAD NADH ATP P2G P3G PEP PYR GLYCERAL3P ADP
clear AMP DHAP GLYC3P GLYCEROL ETOH G1P UTP UDP UDP_GLC TRE T6P PI IMP INO HYP
clear ETOHec GLYCec FRCi