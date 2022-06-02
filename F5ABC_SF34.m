% % Knock-outs to generate
% After the talk with Koen, with the current model, generate in silico 
% mutants: (1) what if the trehalose cycle there was knocked down?, (2) 
% What if glycerol was knocked down, (3) What if both?, (4) What if only 
% one the trehalases (if NTH1, link to dataset FF04).


%% (0) Initial setup.

% startup
clear, %close all
set_paths;
dbstop if error

% assays to run
setup.runSScanelas = 0;
setup.runGSvanHeerden = 0;
setup.plotResultsMode = 0; %30; %20; %10;
setup.conditionsSS = 0;
setup.conditionsGP = 0;

% model structure used
setup.biomass = 0; % biomass
setup.clamp10.Pi = 0;%1; 
setup.clamp10.TRE = 0;%1; 
setup.clamp10.NADX = 0;%1; 
setup.clamp10.AXP = 0;%1; 
setup.clamp10.IXP = 1; 
setup.clamp10.UDP_GLC = 1; % UDP_GLC_exp in AGT1 rate equation
setup.clamp = setup.clamp10;

% select parameters to study:
parsNTH1 = [122, 123];
parsATH1 = [149, 150];
parsAGT1 = [151:154, 158];  % include par158 (Ki for UDPG)
% parsvacuoleT = [155, 156];
parsvacuoleT = [155, 156, 157];
parsGLT = [35:36];
parsGLK = [28, 29, 30, 31, 32, 33, 34]; % 33 is k_T6P
parsPGI = [57:60];
parsPFK = [43:56,87];
parsPGM1 = [83:86];
parsUGP = [144:148];
parsTPS1 = [124:128];
parsTPS2 = [119:121];
parsPFK_sweep = [43:45,47,51,52,87];
parsPFK_selected = [51, 52];
parCombs = [33, parsNTH1, parsATH1, parsAGT1, parsvacuoleT, parsPGM1, parsUGP, parsTPS1, parsTPS2];

% get datasets (others needed)
[legenda] = legendaFull;
loadData_Y3M1;
% get GPFF datasets 
load('TUDdata.mat'); % loads the datasets
    dataset.FF01.time_mets = [0;5;10;15;20;25;30;60;90;120;150;180;220;250;300;350;400];
    dataset.FF01.timeECgluc= [0;5;11;15;20;   30;60;90;    150;180;220;250;300;350;400];
    dataset.FF03.time_mets = [0;5;10;20;30;40;60;90;120;150;200;250;300;400;550;700;800;900;1000;1200;1400;1600;1700;1803];
    dataset.FF04.time_mets = [0;5;10;15;20;30;60;90;120;150;180;220;250;300;350;398];
reorganiseTUDdata; % puts it in order for easily coding in the next section


%% (1) simulations (GP suarez)
% simulation FF01
setup.GPdataset.GP400WT = 1;
setup.GPdataset.GP1800WT = 0;
setup.GPdataset.GP400M = 0;
setup.GPdataset.Fructose  = 0;
setup.GPdataset.Sucrose = 0;

su = 0;
fn = fieldnames(setup.GPdataset);
for i = 1:numel(fn)
    su = su+ setup.GPdataset.(fn{i});
end
if su > 1
    disp('too many datasets used')
end

z = 2; %l/gcdw
if setup.GPdataset.GP400WT        == 1
    InitCond_GP_TUD;
elseif setup.GPdataset.GP1800WT 	== 1
    InitCond_GP_TUD_1800;
elseif setup.GPdataset.GP400M     == 1
    InitCond_GP_TUD_mutant;
elseif setup.GPdataset.Fructose == 1
    FFFruc = readFructose(z);
    dataset.FFFruc = FFFruc; 
    reorganiseFrucData;
    IC = ICFrucData(z, FFFruc);
elseif setup.GPdataset.Sucrose == 1
    FFSuc = readSucrose(z);
    dataset.FFSuc = FFSuc;
    reorganiseSucData;
    InitConditionsSucrose;
end
NumberCycles = 5;
plotflag = 0;


%% updated parameter set from pE2
% % % % % load('pset_pE2res.mat');
% % % % % x1 = x_pE2_start;
% % % % % x2 = x_pE2_end;
%% glycerol_SynthDeg active
setup.glycSynthDeg = 1; % this adds a reaction to the mass balance (ODEs), reaction calculation in rate equations and the starting parameter value
% % % % % x1(159) = -4; % -10;
% % % % % x2(159) = -4; %-10;
%% extension glycogen metabolism
% load parameters
% % % % % load('pset_pE6_x1res.mat','x_pE6_x1_end'); x1 = x_pE6_x1_end;
% % % % % load('pset_pE6_x2res.mat','x_pE6_x2_end'); x2 = x_pE6_x2_end;
%% adding a more detailed glycogen metabolism
setup.glycSynthDeg_separate = 1;
%% recall last case
% % % % % x = x1;
% % % % % x(162) = 10;
setup.ATHinhibitionT6P = 1;
% % % % % load('pset_pE7_xres.mat','x_pE7_start','x_pE7_end');
% % % % % clear x x1 x2
% % % % % x = x_pE7_end;
%% separating ath1_ec and ath1_vac
% % % % % x(163:164) = x(149:150); % km_tre, k_cat
% % % % % x(165) = x(162); % ki_t6p
setup.ATH_separate_EC_VAC = 1;
%% setting the option
setup.glycogenReactionsSink = 1;
setup.dataset = dataset;
%
% % % % % load('pset_pE10_xres.mat','x_pE10_start','x_pE10_end'); x = x_pE10_end;
%% latest setup
% % % % setup.updated_bmf_Cx_ATH1ec = 1;
setup.TREec_brothOut_OFF = 1;
setup.updated_bmf_Cx_ATH1ec = 1;

%% added for enrichment simulations
load('datasetEnrich.mat');
reorganiseEnrichData;

%% (2021 - 08 - 09) right balance around ATH1ec
setup.TRE_recirculation_rightBalances = 1;

%% (2021 09 17) Adjustment glk
% load('x16a.mat','x16a')
% x = x16a;
% x16c_E_start = x_test(1,:);
% x16c_E_final = x_test(2,:);
% load('x16c_E_TPS2.mat', 'x16c_E_start', 'x16c_E_final')
% x = x16c_E_final;
load('pset_Y3M2.mat','x16d_initial','x16d_final')
x = x16d_final;


%% Directly implementing the Csmin idea in out model
% general clamping and setup options
NumberCycles = 5;
setup.clamp_GLCec = 0;
% setup.clamp_GLCec = 1;
setup.changing_Keq_glt = 3; % here the value number 3.
    % changing_Keq_glt = 3; means that CSmin is implemented (FF01 value)
    %                  = 2; manually fixed Keq_glt.  p.GLT.KeqGLT = 1 * 10 .^ setup.Keq_glt_inc;
    %                  = 1; estimated Keq_glt.  p.GLT.KeqGLT = 1 * 10 .^ x(166);
% x(166) = 0; % % <== change x166 if needed here.
setup.csmin = 0.094;
% setup.csmin = 0;

% tic
% %% testing the change in the mass balance
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
legendaMetabolites_addEnrichment;
choosedataset
% 
setup.experiment = 1;


%% (2021 10 01: adding the sinkPYR decrease)
setup.decrease_sinkPYR = 1;
setup.ratio_decrease_sinkPYR = 0.1;


% % % % %% checkup Simulation
% % % % % 
% % % % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % % 
% % % % % T_FF01_1 = T_FF01_2;
% % % % % Y_FF01_1 = Y_FF01_2;
% % % % % V_FF01_1 = V_FF01_2;
% % % % % 
% % % % selResCell{1}.T_FF01 = T_FF01_1;
% % % % selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % selResCell{1}.V_FF01 = V_FF01_1;
% % % % 
% % % % 
% % % % %% Visualization
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % % clear selResCell
% % % % plotMode = 0;


%% (1/8) Recall parameter set Y3M1 and Y3M2
load('pset_Y3M1.mat'); 
x_Y3M1 = x99b;
x_Y3M2 = x; x_Y3M2(35) = x_Y3M2(38); setup.sameKm_Y3M2 = 1; % !!!!!!! (2022 - 01 - 20)
idxs_change = find(x_Y3M1 ~= x_Y3M2(1:length(x_Y3M1)));
%   Columns 1 through 8
% 
%     28    29    30    31    32    33    34    36
% 
%   Columns 9 through 16
% 
%     37    38    83    84    85    86   109   119
% 
%   Columns 17 through 24
% 
%    120   121   122   123   124   125   126   127
% 
%   Columns 25 through 32
% 
%    128   129   144   145   146   147   148   149
% 
%   Columns 33 through 35
% 
%    150   151   152

% Changes which are:
% % % % parsGLT = [35:36,38];
% % % % parsGLK = [28, 29, 30, 31, 32, 33, 34]; % 33 is k_T6P
% % % % parsPGM1 = [83:86];
% % % % parsUGP = [144:148];
% % % % parsTPS1 = [124:128];
% % % % parsTPS2 = [119:121];
% % % % parsNTH1 = [122, 123];
% % % % parsATH1 = [149, 150];
% % % % parsAGT1 = [151:154, 158];  % include par158 (Ki for UDPG)

% % other parameters
% parsPGI = [57:60];
% parsPFK = [43:56,87];
% %
x_Y3M2(109) = x_Y3M1(109); % ATPase changed here directly
x_Y3M2(129) = x_Y3M1(129) + 1; % ATPase changed here directly
x_Y3M1 = [x_Y3M1, zeros(1,13)];
% 
% x2 = x; x2(109) = x_Y3M1(109); x2(129) = x_Y3M1(129)+1;

%% (2/8) Simulate both
% % % % clear selResCell
% % % % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % [T_FF01_2,Y_FF01_2,V_FF01_2] = simulate_FF_enrichment(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % [T_FF01_2,Y_FF01_2,V_FF01_2] = simulate_FF_enrichment(x_Y3M1,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % [T_FF01_3,Y_FF01_3,V_FF01_3] = simulate_FF_enrichment(x_Y3M2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % 
% % % % % 
% % % % selResCell{1}.T_FF01 = T_FF01_1;
% % % % selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % selResCell{1}.V_FF01 = V_FF01_1;
% % % % % 
% % % % selResCell{2}.T_FF01 = T_FF01_2;
% % % % selResCell{2}.Y_FF01 = Y_FF01_2;
% % % % selResCell{2}.V_FF01 = V_FF01_2;
% % % % % 
% % % % selResCell{3}.T_FF01 = T_FF01_3;
% % % % selResCell{3}.Y_FF01 = Y_FF01_3;
% % % % selResCell{3}.V_FF01 = V_FF01_3;
% % % % 
% % % % % %% Visualization
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;

%% (2/8)b. Pre-screen for expected changes
% actual changes
parsGLT = [35:36,38];
parsGLK = [28, 29, 30, 31, 32, 33, 34]; % 33 is k_T6P
parsPGM1 = [83:86];
parsUGP = [144:148];
parsTPS1 = [124:128];
parsTPS2 = [119:121];
parsNTH1 = [122, 123];
parsATH1 = [149, 150];
parsAGT1 = [151:154, 158];  % include par158 (Ki for UDPG)
    parsvacuoleT = [155, 156, 157];
    parsGLYC = 159:161;
    parsath1_T6P = 162;
    parsATH1ec = 163:165;
%     parsAGT1_T6P = 166;
    parsLast = [parsvacuoleT, parsGLYC, parsath1_T6P, parsATH1ec];
parsTRE = [parsPGM1, parsUGP, parsTPS1, parsTPS2, parsNTH1, parsATH1, parsAGT1];
parsTRE2 = [parsTRE, parsLast];
% other parameters
parsPGI = [57:60];
parsPFK = [43:56,87];
parsALD = [11:15];
parsGAPDH = [21:27];

% pars...

x_Y3M1_parsGLT = x_Y3M1; x_Y3M1_parsGLT(parsGLT) = x_Y3M2(parsGLT);
    x_Y3M1_parsGLT(129) = x_Y3M2(129);
    x_Y3M1_parsGLT(37) = x_Y3M2(37);
x_Y3M1_parsGLK = x_Y3M1; x_Y3M1_parsGLK(parsGLK) = x_Y3M2(parsGLK);
    x_Y3M1_parsGLK(129) = x_Y3M2(129);
    x_Y3M1_parsGLK(37) = x_Y3M2(37);
x_Y3M1_parsTRE = x_Y3M1; x_Y3M1_parsTRE(parsTRE2) = x_Y3M2(parsTRE2);
    x_Y3M1_parsTRE(129) = x_Y3M2(129);
    x_Y3M1_parsTRE(37) = x_Y3M2(37);
x_Y3M1_parsGLTGLKTRE = x_Y3M1; x_Y3M1_parsGLTGLKTRE([parsGLT,parsGLK,parsTRE2]) = x_Y3M2([parsGLT,parsGLK,parsTRE2]);
    x_Y3M1_parsGLTGLKTRE(129) = x_Y3M2(129);
    x_Y3M1_parsGLTGLKTRE(37) = x_Y3M2(37);
x_Y3M1_129 = x_Y3M1;
    x_Y3M1_129(129) = x_Y3M2(129);
    x_Y3M1_129(37) = x_Y3M2(37);
    
% % % % % % % % % % simulate
% % % % % % % % % clear selResCell
% % % % % % % % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x_Y3M2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % % % % [T_FF01_2,Y_FF01_2,V_FF01_2] = simulate_FF_enrichment(x_Y3M1_parsGLT,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % % % % [T_FF01_3,Y_FF01_3,V_FF01_3] = simulate_FF_enrichment(x_Y3M1_parsGLK,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % [T_FF01_4,Y_FF01_4,V_FF01_4] = simulate_FF_enrichment(x_Y3M1_parsTRE,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % % % % [T_FF01_5,Y_FF01_5,V_FF01_5] = simulate_FF_enrichment(x_Y3M1_parsGLTGLKTRE,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % % % % [T_FF01_6,Y_FF01_6,V_FF01_6] = simulate_FF_enrichment(x_Y3M1_129,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % % % % %%
% % % % % % % % clear selResCell
% % % % % % % % % x_Y3M1_parsGLT = x_Y3M1; x_Y3M1_parsGLT(parsGLT) = x_Y3M2(parsGLT);
% % % % % % % % % Pretty off.
% % % % % % % % % x_Y3M1_parsGLK = x_Y3M1; x_Y3M1_parsGLK(parsGLK) = x_Y3M2(parsGLK);
% % % % % % % % % Pretty off
% % % % % % % % % x_Y3M1_parsTRE = x_Y3M1; x_Y3M1_parsTRE(parsTRE2) = x_Y3M2(parsTRE2);
% % % % % % % % % getting there
% % % % % % % % % % x_Y3M1_parsGLTGLKTRE = x_Y3M1; x_Y3M1_parsGLTGLKTRE([parsGLT,parsGLK,parsTRE2]) = x_Y3M2([parsGLT,parsGLK,parsTRE2]);
% % % % % % % % % Same as initial.
% % % % % % % % % % x_Y3M1_129 = x_Y3M1;
% % % % % % % % % Pretty off, as expected.
% % % % % 
% % % % selResCell{1}.T_FF01 = T_FF01_4;
% % % % selResCell{1}.Y_FF01 = Y_FF01_4;
% % % % selResCell{1}.V_FF01 = V_FF01_4;
% % % % % % % % % 
% % % % % % % % selResCell{1}.T_FF01 = T_FF01_1;
% % % % % % % % selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % % % % % selResCell{1}.V_FF01 = V_FF01_1;
% % % % % % % % % 
% % % % % % % % selResCell{2}.T_FF01 = T_FF01_6;
% % % % % % % % selResCell{2}.Y_FF01 = Y_FF01_6;
% % % % % % % % selResCell{2}.V_FF01 = V_FF01_6;
% % % % % % % % % 
% % % % % % % % selResCell{3}.T_FF01 = T_FF01_4;
% % % % % % % % selResCell{3}.Y_FF01 = Y_FF01_4;
% % % % % % % % selResCell{3}.V_FF01 = V_FF01_4;
% % % % % % % % % % 
% % % % % % % % % selResCell{2}.T_FF01 = T_FF01_2;
% % % % % % % % % selResCell{2}.Y_FF01 = Y_FF01_2;
% % % % % % % % % selResCell{2}.V_FF01 = V_FF01_2;
% % % % % % % % % % 
% % % % % % % % % selResCell{3}.T_FF01 = T_FF01_3;
% % % % % % % % % selResCell{3}.Y_FF01 = Y_FF01_3;
% % % % % % % % % selResCell{3}.V_FF01 = V_FF01_3;
% % % % % % % % % % 
% % % % % % % % % selResCell{4}.T_FF01 = T_FF01_4;
% % % % % % % % % selResCell{4}.Y_FF01 = Y_FF01_4;
% % % % % % % % % selResCell{4}.V_FF01 = V_FF01_4;
% % % % % % % % % % 
% % % % % % % % % selResCell{5}.T_FF01 = T_FF01_5;
% % % % % % % % % selResCell{5}.Y_FF01 = Y_FF01_5;
% % % % % % % % % selResCell{5}.V_FF01 = V_FF01_5;
% % % % % % % % % % 
% % % % % % % % % selResCell{6}.T_FF01 = T_FF01_6;
% % % % % % % % % selResCell{6}.Y_FF01 = Y_FF01_6;
% % % % % % % % % selResCell{6}.V_FF01 = V_FF01_6;
% % % % % % % % 
% % % % % %% Visualization
% % % % plotMode = 0;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;
% % % % 
% % % % 
%%
[T_FF01_ref,Y_FF01_ref,V_FF01_ref] = simulate_FF_enrichment(x_Y3M2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% 
selResCell_ref{1}.T_FF01 = T_FF01_ref;
selResCell_ref{1}.Y_FF01 = Y_FF01_ref;
selResCell_ref{1}.V_FF01 = V_FF01_ref;
% 
% setup to estimate from the reference simulation
% setup.simData = selResCell_ref;
setup.simData = selResCell_ref{1};
setup.pset_Y3M1 = x_Y3M1;


%%
clear selResCell
selResCell = selResCell_ref;
% 
plotMode = 0;
referencePlotSimulations_enrichment
plotMode = 0;


%% (2/8)b. Pre-screen for expected changes
% 
x_Y3M1_parsTRE_GLT = x_Y3M1; x_Y3M1_parsTRE_GLT([parsGLT, parsTRE2]) = x_Y3M2([parsGLT, parsTRE2]);
    x_Y3M1_parsTRE(129) = x_Y3M2(129);
    x_Y3M1_parsTRE(37) = x_Y3M2(37);
% 
x_Y3M1_parsTRE_GLK = x_Y3M1; x_Y3M1_parsTRE_GLK([parsGLK, parsTRE2]) = x_Y3M2([parsGLK, parsTRE2]);
    x_Y3M1_parsTRE(129) = x_Y3M2(129);
    x_Y3M1_parsTRE(37) = x_Y3M2(37);

% % % % % % simulate
% % % % % clear selResCell
% % % % [T_FF01_11,Y_FF01_11,V_FF01_11] = simulate_FF_enrichment(x_Y3M2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % [T_FF01_12,Y_FF01_12,V_FF01_12] = simulate_FF_enrichment(x_Y3M1_parsTRE_GLT,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % [T_FF01_13,Y_FF01_13,V_FF01_13] = simulate_FF_enrichment(x_Y3M1_parsTRE_GLK,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % 
% % % % % %%
% % % % clear selResCell
% % % % % 
% % % % selResCell{1}.T_FF01 = T_FF01_11;
% % % % selResCell{1}.Y_FF01 = Y_FF01_11;
% % % % selResCell{1}.V_FF01 = V_FF01_11;
% % % % % 
% % % % selResCell{2}.T_FF01 = T_FF01_12;
% % % % selResCell{2}.Y_FF01 = Y_FF01_12;
% % % % selResCell{2}.V_FF01 = V_FF01_12;
% % % % % 
% % % % selResCell{3}.T_FF01 = T_FF01_13;
% % % % selResCell{3}.Y_FF01 = Y_FF01_13;
% % % % selResCell{3}.V_FF01 = V_FF01_13;
% % % % 
% % % % 
% % % % % %% Visualization
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;


%% (3/8) Prepare the multiple combinations of parameters
% blank and constant setup
legendaMetabolites_addEnrichment;
choosedataset
setup.ExpData = ExpData;
% 
blankWeight = zeros(1,89); % 85+3 for glycerol, 89 for glucose enrichment
lambdalist = 0;
setup.parEst.lambda = lambdalist(1); lam = setup.parEst.lambda;
% % % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx});
% % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
% %     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% options = optimoptions('lsqnonlin','Display','iter',...
%     'OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
%     'FiniteDifferenceStepSize',0.2);
options = optimoptions('lsqnonlin','Display','iter',...
    'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4,...
    'DiffMinChange',0.1);
%     'OutputFcn',{@saveIterationsMain},...
%     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% setup.parEst.costfun = 10015; 

% creating enzyme combinations
enzNum = 6;
% 1 -> GLT
% 2 -> GLK
% 3 -> PGI
% 4 -> PFK
% 5 -> ALD
% 6 -> GAPDH
nchoosek_cell = cell(enzNum,1);
for i = 1:enzNum
    nchoosek_cell{i} = nchoosek(1:enzNum,i);
end
% concatenate
if exist('z','var')
    clear z
end
for i = 1:(enzNum-1)
    if exist('z','var')
    else
        z = nchoosek_cell{i}';
    end
    B = nchoosek_cell{i+1}';
    sA = size(z);
    sB = size(B);
    C = max(sA(1),sB(1));
    z = [[z;zeros(abs([C 0]-sA))],[B;zeros(abs([C,0]-sB))]];
end
z = z';

% number of tests + initial parameter set
% z = z([1:3,5:6,8,11,15],:);
[nTests,nZ] = size(z);
xStart = x_Y3M1; 
    xStart([37,129]) = x_Y3M2([37,129]);
    xStart(parsTRE2) = x_Y3M2(parsTRE2);
% x_Y3M1_parsTRE
% % names_cell
names_cell = cell(nTests,1);
for i = 1:nTests
    names_cell{i} = sprintf('FF_mFig_7C_n%d.mat',i);
end

% % parComb_cell
parComb_cell = cell(nTests,1);
for i = 1:nTests
%     tempPars = parsTRE2;
    tempPars = [];
    if any(z(i,:) == 1), tempPars = [tempPars, parsGLT]; end
    if any(z(i,:) == 2), tempPars = [tempPars, parsGLK]; end
    if any(z(i,:) == 3), tempPars = [tempPars, parsPGI]; end
    if any(z(i,:) == 4), tempPars = [tempPars, parsPFK]; end
    if any(z(i,:) == 5), tempPars = [tempPars, parsALD]; end
    if any(z(i,:) == 6), tempPars = [tempPars, parsGAPDH]; end
    parComb_cell{i} = tempPars;
end

% % cost_function
% improving 50% for reaction rates (all of UG) and 50% for metabolites (UG and TRE-cycle)
w_glc_ic = 6;
w_glc_ec = 36; % concentrations
w_g1p = 21; % concentrations
w_udpg = 24; % concentrations
w_t6p_tre = [26 25]; % concentrations
w_treRates_glk_pgi = [55 59 57 58 56 40 41]; % rates
w_glt = 39; % rates
w_gly = [86 87 88]; % concentrations
w_tre_ic = 25;
w_tre_ec = 37; % concentrations
w_ath1 = [83 84];
    w_ath1_ec = 83;
    w_ath1_vac = 84; 
w_agt1 = 85;
w_glc_enrich = 89;
% concentrations and rates
wConcentrations = [w_glc_ec, w_g1p, w_udpg, w_t6p_tre, 86, w_tre_ec];
wRates = [w_treRates_glk_pgi, w_glt];
% adding to the cost function matrix

% % warray_cell
warray_cell = cell(nTests,1);
for i = 1:nTests % base,
    warray_cell{i} = blankWeight;
%     warray_cell{i}(wConcentrations) = 0.5 * ones / length(wConcentrations);
%     warray_cell{i}(wRates) = 0.5 * ones / length(wRates);
    warray_cell{i}([62 65]) = ones;
end


% % % % %% (4/8) Run them
% % % % % parSet
% % % % % parSet1 = 1:2:63;
% % % % % parSet2 = 2:2:62;
% % % % parSet1 = 1:40;
% % % % parSet2 = 41:63;
% % % % % for o = selectedRuns_idx
% % % % cluster = parcluster('local');
% % % % % pool = parpool(cluster,4);
% % % % pool = parpool(cluster,16);
% % % % parfor o = 1:nTests
% % % % % parfor o = parSet1
% % % % % parfor o = parSet2
% % % % % for o = 1:nTests
% % % % 
% % % %     % disp
% % % %     disp(o);
% % % % %     fprintf('o = %f.\n', o) %disp(o)
% % % % 
% % % %     % find in the precreated arrays
% % % %     selPars = parComb_cell{o};
% % % %     warray = warray_cell{o}; 
% % % %     saveName = names_cell{o};
% % % % %     setup.selPars = selPars;
% % % %     
% % % %     % core options:
% % % %     plength = length(selPars);
% % % %     lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
% % % %     ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;
% % % %     x_temp = xStart(selPars);
% % % % 
% % % % %     % %% run check
% % % % %     [error]=costfunSystemY3M1_FF_pRun_enrichment(x_temp,canelas_SS,setup,xStart,data,dataset,NumberCycles,IC,selPars,warray);
% % % % %     % errorAnalysis_Y3M1;
% % % % %     % %%
% % % % % 
% % % % %     %
% % % % %     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(xStart,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % %     plotMode = 2; % multiple simulations
% % % % %     % clear selResCell
% % % % %     selResCell{1}.T_FF01 = T_FF01_1;
% % % % %     selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % %     selResCell{1}.V_FF01 = V_FF01_1;
% % % % %     % 
% % % % %     referencePlotSimulations_enrichment
% % % % %     plotMode = 0;
% % % % %     %
% % % %         
% % % %         % parameter estimation
% % % % %         NumberCycles = 20;
% % % %         sprintf('run num.%d tested',o)
% % % %         tic
% % % % %         [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun,x_temp,lb,ub,options,canelas_SS,setup,x_iter,data,dataset,NumberCycles,IC,selPars,warray);
% % % %         [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun_enrichment,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,NumberCycles,IC,selPars,warray);
% % % %         t = toc; 
% % % %         disp(xres);
% % % %         
% % % %         % simulating output
% % % %         xRound = xStart;
% % % %         xRound(selPars) = xres;
% % % %         
% % % %         [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xRound,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %         simRes1_T_FF01 = T_FF01;
% % % %         simRes1_Y_FF01 = Y_FF01;
% % % %         simRes1_V_FF01 = V_FF01;
% % % %         
% % % %         %         plotMode = 2; % multiple simulations
% % % % %         % clear selResCell
% % % % %         selResCell{1}.T_FF01 = T_FF01_1;
% % % % %         selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % %         selResCell{1}.V_FF01 = V_FF01_1;
% % % %         
% % % % 
% % % %         % saving
% % % % %         parsave_Y3M2_cluster(saveName, xres, resnorm, residual, exitflag, warray, t, selPars);
% % % %         parsave_Y3M2_cluster2(saveName, xres, resnorm, residual, ...
% % % %             exitflag, warray, t, selPars, ...
% % % %             simRes1_T_FF01, simRes1_Y_FF01, simRes1_V_FF01);
% % % % %         save(saveName, 'xres', 'resnorm', 'residual', 'exitflag', 'warray', 't', 'selPars');
% % % %     
% % % % end
% % % % delete(pool)
% % % % quit


% %% (5/8) Recall results, visualize and select (Option 2)
% % 
% xAll1 = xStart;
% %
% namesHits1 = cell(1,1); namesHits1{1} = 'initial';
% % 
% idxs = 0;
% %
% simRes1 = selResCell_ref;
% errorRes1 = cell(1,1); errorRes1{1} = zeros(409,1);
% %
% for i = 1:nTests
%     % pC1_wA1
% %     temp_idx = i;
% %     loadName = sprintf(names_cell{temp_idx}, i);
%     loadName = names_cell{i};
%     if exist(loadName,'file') == 2 % if exist
%         load(loadName);
%         selPars = parComb_cell{i};
%         x3 = xStart;
%         x3(selPars) = xres;
%         xAll1 = [xAll1; x3];
%         namesHits1 = [namesHits1; loadName];
%         idxs = [idxs, i];
%         %
%         temp_simRes.T_FF01 = simRes1_T_FF01;
%         temp_simRes.Y_FF01 = simRes1_Y_FF01;
%         temp_simRes.V_FF01 = simRes1_V_FF01;
%         simRes1 = [simRes1, temp_simRes];
%         tempErrCell = cell(1,1); tempErrCell{1} = residual;
%         errorRes1 = [errorRes1; tempErrCell];
%         %
%     end
% end
% % recall
% simRes1_mF7_sims1 = simRes1;
% simRes1_mF7_errorRes1 = errorRes1;
% simRes1_16e_xAll1 = xAll1;
% simRes1_16e_namesHits1 = namesHits1;
% % save
% save('tempRes2_mF7CnewStart.mat', ...
%     'simRes1_mF7_sims1', 'simRes1_mF7_errorRes1', ...
%     'simRes1_16e_xAll1', 'simRes1_16e_namesHits1');


%%
% % % First save
% % load('tempRes2_mF7.mat', ...
% %     'simRes1_mF7_sims1', 'simRes1_mF7_errorRes1', ...
% %     'simRes1_16e_xAll1', 'simRes1_16e_namesHits1');
% % Second save
% load('tempRes2_mF7_b.mat', ...
%     'simRes1_mF7_sims1', 'simRes1_mF7_errorRes1', ...
%     'simRes1_16e_xAll1', 'simRes1_16e_namesHits1');
% Second save
load('tempRes2_mF7CnewStart.mat', ...
    'simRes1_mF7_sims1', 'simRes1_mF7_errorRes1', ...
    'simRes1_16e_xAll1', 'simRes1_16e_namesHits1');


%% Visualizing results: simulation
% selResCell = simRes1_16e_TRE([1,11:19]);
clear selResCell
selResCell = [selResCell_ref];
selResCell = simRes1_mF7_sims1(1:7);
% selResCell = simRes1_mF7_sims1(1:8);
% selResCell = [selResCell_ref, simRes1_mF7_sims1(2:4)];
% selResCell = simRes1_mF7_sims1(1:3);
% selResCell = simRes1_mF7_sims1(8:12);
% selResCell = simRes1_mF7_sims1;

plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;


%% calculate the error
% calculate error 2, based on error rates glycolysis + trehalsoe cycle

% recall idxs
% refSim
time2compare = [0:2:30,40:10:160,180:20:400]; % time2compare
A = repmat(selResCell_ref{1}.T_FF01,[1 length(time2compare')]);
[~,closestIndex] = min(abs(A-time2compare));
% N(closestIndex) % find closest value
vGLTref = selResCell_ref{1}.V_FF01(closestIndex,1);
vGLKref = selResCell_ref{1}.V_FF01(closestIndex,2);
vPGIref = selResCell_ref{1}.V_FF01(closestIndex,3);
vPFKref = selResCell_ref{1}.V_FF01(closestIndex,4);
vALDref = selResCell_ref{1}.V_FF01(closestIndex,5);
vPGM1ref = selResCell_ref{1}.V_FF01(closestIndex,17);
vUGPref = selResCell_ref{1}.V_FF01(closestIndex,18);
vTPS1ref = selResCell_ref{1}.V_FF01(closestIndex,21);
vTPS2ref = selResCell_ref{1}.V_FF01(closestIndex,19);
vNTH1ref = selResCell_ref{1}.V_FF01(closestIndex,20); % fluxes ref
    % % 
    % figure, 
    % plot(selResCell_ref{1}.T_FF01(closestIndex),vTPS2ref,'.-')
% blaks for the other sims
errorTre = zeros(length(simRes1_16e_namesHits1) - 1, 1);    % make blank cell error trehalose 
errorGly = zeros(length(simRes1_16e_namesHits1) - 1, 1);    % make blank cell error glycolysis 
errorTotal = zeros(length(simRes1_16e_namesHits1) - 1, 1);    % make blank cell error total 
idxs2 = zeros(length(simRes1_16e_namesHits1) - 1, 1);        % make blank idxs2 
numEnzymes = zeros(length(simRes1_16e_namesHits1) - 1, 1);    % make blank num enzumes 
enzNames = {'HXT','GLK','PGI','PFK','ALD','TDH'};
cellNamesCombs = cell(length(simRes1_16e_namesHits1) - 1, 1);    % make blank cell names combinations 
% otherSims
for i = 2:length(simRes1_16e_namesHits1) % other simulations (skip initial point)
% for i = 2:3
    A = repmat(simRes1_mF7_sims1{i}.T_FF01,[1 length(time2compare')]);
    [~,closestIndex] = min(abs(A-time2compare)); % find closest value
    % fluxes exp
    vGLTsim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,1);
    vGLKsim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,2);
    vPGIsim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,3);
    vPFKsim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,4);
    vALDsim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,5);
    vPGM1sim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,17);
    vUGPsim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,18);
    vTPS1sim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,21);
    vTPS2sim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,19);
    vNTH1sim = simRes1_mF7_sims1{i}.V_FF01(closestIndex,20);
%     % fluxes error (method 2: denominator = DATAsim)
%     vGLTerror = sum(abs((vGLTsim - vGLTref) ./ vGLTsim));
%     vGLKerror = sum(abs((vGLKsim - vGLKref) ./ vGLKsim));
%     vPGIerror = sum(abs((vPGIsim - vPGIref) ./ vPGIsim));
%     vPFKerror = sum(abs((vPFKsim - vPFKref) ./ vPFKsim));
%     vALDerror = sum(abs((vALDsim - vALDref) ./ vALDsim));
%     vPGM1error = sum(abs((vPGM1sim - vPGM1ref) ./ vPGM1sim));
%     vUGPerror = sum(abs((vUGPsim - vUGPref) ./ vUGPsim));
%     vTPS1error = sum(abs((vTPS1sim - vTPS1ref) ./ vTPS1sim));
%     vTPS2error = sum(abs((vTPS2sim - vTPS2ref) ./ vTPS2sim));
%     vNTH1error = sum(abs((vNTH1sim - vNTH1ref) ./ vNTH1sim));
    % fluxes error (method 2: denominator = max and single value)
    vGLTerror = sum(abs((vGLTsim - vGLTref) ./ max(vGLTref)));
    vGLKerror = sum(abs((vGLKsim - vGLKref) ./ max(vGLKref)));
    vPGIerror = sum(abs((vPGIsim - vPGIref) ./ max(vPGIref)));
    vPFKerror = sum(abs((vPFKsim - vPFKref) ./ max(vPFKref)));
    vALDerror = sum(abs((vALDsim - vALDref) ./ max(vALDref)));
    vPGM1error = sum(abs((vPGM1sim - vPGM1ref) ./ max(vPGM1ref)));
    vUGPerror = sum(abs((vUGPsim - vUGPref) ./ max(vUGPref)));
    vTPS1error = sum(abs((vTPS1sim - vTPS1ref) ./ max(vTPS1ref)));
    vTPS2error = sum(abs((vTPS2sim - vTPS2ref) ./ max(vTPS2ref)));
    vNTH1error = sum(abs((vNTH1sim - vNTH1ref) ./ max(vNTH1ref)));
%         %
%         figure
%         subplot(1,2,1)
%         plot(simRes1_mF7_sims1{i}.T_FF01(closestIndex), vTPS2ref, '.-')
%         hold on
%         plot(simRes1_mF7_sims1{i}.T_FF01(closestIndex), vTPS2sim, '.-')
%         legend('ref','sim')
%         subplot(1,2,2)
%         plot(simRes1_mF7_sims1{i}.T_FF01(closestIndex), (vTPS2sim - vTPS2ref) ./ 1, '.-')
%         legend('error')
    % get the errors
%     errorGly(i-1) = vGLTerror + vGLKerror + vPGIerror + vPFKerror + vALDerror;    % to calculate the trehalose errors (weighted)
    errorGly(i-1) = vGLTerror;    % to calculate the trehalose errors (weighted)
%     errorTre(i-1) = vPGM1error + vUGPerror + vTPS1error + vTPS2error + vNTH1error;    % to calculate the glycolysis errors (weighted)
    errorTre(i-1) = vTPS2error;    % to calculate the glycolysis errors (weighted)
    errorTotal(i-1) = errorTre(i-1) + errorGly(i-1);    % to calculate the total errors (weighted)
    idxs2(i-1) = str2double(erase(erase(simRes1_16e_namesHits1{i},".mat"),"FF_mFig_7C_n"));    % update idxs2 at each loop (detecit the number in simRes1_16_namesHits1
    numEnzymes(i-1) = sum(abs(z(idxs2(i-1),:) ~= 0));    % calculate num error (sum is_num)
    tempNamesCombs = [];
    for j = 1:enzNum
        if ismember(j,z(idxs2(i-1),:))
            tempNamesCombs = [tempNamesCombs, enzNames{j}, ', '];
        end
    end
    cellNamesCombs{i-1} = tempNamesCombs;   % recall name enzyme combinations (i-1)
end
% 


%%
figure,
subplot(1,3,1), plot(1:length(idxs2), errorGly), title('errorGly')
subplot(1,3,2), plot(1:length(idxs2), errorTre), title('errorTre')
subplot(1,3,3), plot(1:length(idxs2), errorTotal), title('errorTotal')


%% visual (1/4) error vs numPars, barplots
% plot per option, as made before

% barplot, orientation vertical
% at first, plot on top errorGLYCOLYSIS and errorTREcycle
% add the plot of the number of enzumes
tempName = categorical(cellNamesCombs);
tempName = reordercats(tempName,cellNamesCombs);

% clf(101)
fig101_h = figure(101);
barh(tempName, [errorGly, errorTre],...
    'stacked'), 
legend('errorGly','errorTre')
ax = gca;
ax.FontSize = 8; 


%% Making the figure plot here
if exist('figh','var')
    clf(1001)
end
figh = figure(1001);
figh.Position = [962 42 958 954];
% 
% sph = subplot(3,2,[1 3]);
sph = subplot(2,2,[1 3]);
bplot = barh(tempName, errorTotal, 'stacked','facecolor', 'flat');
% bplot.BarWidth = 0.5;
bplot.BarWidth = 0.65;
% bplot.FaceColor = 'black';
% bplot.EdgeColor = 'none';
bplot.EdgeColor = 'black';
sph.XLim = [0 30];
ax = gca;
% ax.FontSize = 6; 
set(gca,'YTickLabel',[]);

% 
hold on
% 
thesis_colopalette
%

% % % % % % % % % % % % % % % % % % 
% Get the data for all the bars that were plotted
x_1 = get(bplot,'XData');
y_1 = get(bplot,'YData');

ygap = 1;  % Specify vertical gap between the bar and label
% ylimits = get(gca,'YLim');
xlimits = get(gca,'XLim');
set(gca,'XLim',[xlimits(1),xlimits(2)+0.2*max(y_1)]); % Increase y limit for labels

% Create labels to place over bars
labels = cellNamesCombs;
for i = 1:length(x_1) % Loop over each bar 
    xpos = x_1(i);        % Set x position for the text label
    ypos = y_1(i) + ygap; % Set y position, including gap
    htext = text(-0.2,xpos,labels{i});          % Add text label
    htext.Position = htext.Position + [0 -.7 0];
    set(htext,'VerticalAlignment','bottom',...  % Adjust properties
        'HorizontalAlignment','right',...
        'FontSize',8)
    % 
    if contains(htext.String, 'HXT, GLK, ')
        set(htext, 'color', cpal_broad.blue)
        bplot.CData(i,:) = cpal_broad.blue;
    else
        bplot.CData(i,:) = [0 0 0];
    end
end
% re adjustment axis
sph.Position = sph.Position + [0.05 0 -0.05 0];

% 
yyaxis right
line([5 5],[0 1],'Color','none','LineStyle','--')
hold on,
% 
% % % % line([5 5],[0 0.095],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([5 10],[0.095 0.095],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([10 10],[0.095 0.3325],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([10 15],[0.3325 0.3325],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([15 15],[0.3325 0.65],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([15 20],[0.65 0.65],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([20 20],[0.65 0.89],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([20 25],[0.89 0.89],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([25 25],[0.89 0.985],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([25 30],[0.985 0.985],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)
% % % % line([30 30],[0.985 1],'Color',cpal_broad.vermillion,'LineStyle','-','LineWidth',2)



%
annotation('arrow',[0.05 0.115],[0.193 0.193],'color',cpal_broad.blue)
% 
set(gca,'YTickLabel',[],'YColor', 'k');%,'color','k');
% 
% % % % text(4.5,-0.015,'1','color',cpal_broad.vermillion)
% % % % text(9.5,-0.015,'2','color',cpal_broad.vermillion)
% % % % text(14.5,-0.015,'3','color',cpal_broad.vermillion)
% % % % text(19.5,-0.015,'4','color',cpal_broad.vermillion)
% % % % text(24.5,-0.015,'5','color',cpal_broad.vermillion)
% % % % text(29.5,-0.015,'6','color',cpal_broad.vermillion)

% 
% % % % temp = text(17.5,-0.05,'Number of enzymes','color',cpal_broad.vermillion,...
% % % %     'FontSize',12,'HorizontalAlignment','center');
% % % % temp2 = text(17.5,1.05,'Model error','color','black',...
% % % %     'FontSize',12,'HorizontalAlignment','center');
% % % % set(gca,'xaxisLocation','top')
temp2 = text(17.5,-0.05,'Model error','color','black',...
    'FontSize',12,'HorizontalAlignment','center');

%% Regularization plot
% To consider changes from here.




%% subplot(2,2,2) -> lambdalist, with pareto front inside
% recall data pars
    % -> which is just a backbone, estimated parameters below
load('pset_pE10_xres.mat','x_pE10_start','x_pE10_end'); %x10 = x_pE10_end;
xAll2 = x_pE10_end;
for warray_16a = 1
    %
    namesHits2 = cell(1,1); namesHits2{1} = 'initial';
    % 
    names_cell2{1} = 'FF_pE16a_glk_wA1.mat';
    names_cell2{2} = 'FF_pE16a_glk_wA2.mat';
    names_cell2{3} = 'FF_pE16a_glk_wA3.mat';
    names_cell2{4} = 'FF_pE16a_glk_wA4.mat';
    names_cell2{5} = 'FF_pE16a_glk_wA5.mat';
    names_cell2{6} = 'FF_pE16a_glk_wA6.mat';
    names_cell2{7} = 'FF_pE16a_glk_wA7.mat';
    names_cell2{8} = 'FF_pE16a_glk_wA8.mat';
    %
    names_cell2{9} = 'FF_pE16a_glk_wA9.mat';
    names_cell2{10} = 'FF_pE16a_glk_wA10.mat';
    names_cell2{11} = 'FF_pE16a_glk_wA11.mat';
    names_cell2{12} = 'FF_pE16a_glk_wA12.mat';
    names_cell2{13} = 'FF_pE16a_glk_wA13.mat';
    names_cell2{14} = 'FF_pE16a_glk_wA14.mat';
    names_cell2{15} = 'FF_pE16a_glk_wA15.mat';
    names_cell2{16} = 'FF_pE16a_glk_wA16.mat';
    %
    parsHXK = 28:34; % hxk
    pC16a = parsHXK;
    % pC16b = [parsGLT, parsHXK, parsTreCycle];
    %
    parComb_cell2 = cell(16,1);
    % 
    parComb_cell2{1} = pC16a;
    parComb_cell2{2} = pC16a;
    parComb_cell2{3} = pC16a;
    parComb_cell2{4} = pC16a;
    parComb_cell2{5} = pC16a;
    parComb_cell2{6} = pC16a;
    parComb_cell2{7} = pC16a;
    parComb_cell2{8} = pC16a;
    parComb_cell2{9} = pC16a;
    parComb_cell2{10} = pC16a;
    parComb_cell2{11} = pC16a;
    parComb_cell2{12} = pC16a;
    parComb_cell2{13} = pC16a;
    parComb_cell2{14} = pC16a;
    parComb_cell2{15} = pC16a;
    parComb_cell2{16} = pC16a;
    %
    for i = 1:16
        % pC1_wA1
        temp_idx = i;
        loadName = sprintf(names_cell2{temp_idx}, i);
        if exist(loadName,'file') == 2 % if exist
            load(loadName);
            selPars = parComb_cell2{temp_idx};
            x3 = x_pE10_end;
            x3(selPars) = xres;
            xAll2 = [xAll2; x3];
            namesHits2 = [namesHits2; loadName];
        end
    end
end
for warray_16a2 = 1
    % names_cell
    names_cell = cell(16,1);
    names_cell{1} = 'FF_pE16a2_glk_wA1.mat';
    names_cell{2} = 'FF_pE16a2_glk_wA2.mat';
    names_cell{3} = 'FF_pE16a2_glk_wA3.mat';
    names_cell{4} = 'FF_pE16a2_glk_wA4.mat';
    names_cell{5} = 'FF_pE16a2_glk_wA5.mat';
    names_cell{6} = 'FF_pE16a2_glk_wA6.mat';
    names_cell{7} = 'FF_pE16a2_glk_wA7.mat';
    names_cell{8} = 'FF_pE16a2_glk_wA8.mat';
    names_cell{9} = 'FF_pE16a2_glk_wA9.mat';
    names_cell{10} = 'FF_pE16a2_glk_wA10.mat';
    names_cell{11} = 'FF_pE16a2_glk_wA11.mat';
    names_cell{12} = 'FF_pE16a2_glk_wA12.mat';
    names_cell{13} = 'FF_pE16a2_glk_wA13.mat';
    names_cell{14} = 'FF_pE16a2_glk_wA14.mat';
    names_cell{15} = 'FF_pE16a2_glk_wA15.mat';
    names_cell{16} = 'FF_pE16a2_glk_wA16.mat';
    
    % % second   
    xAll1 = x;
    %
    namesHits1 = cell(1,1); namesHits1{1} = 'initial';
    %
    for i = 1:16
        % pC1_wA1
        temp_idx = i;
        loadName = sprintf(names_cell{temp_idx}, i);
        if exist(loadName,'file') == 2 % if exist
            load(loadName);
            selPars = parComb_cell2{temp_idx};
            x3 = x;
            x3(selPars) = xres;
            xAll1 = [xAll1; x3];
            namesHits1 = [namesHits1; loadName];
        end
    end
end


% %% recall data sims
% recall first round
load('tempRes_16a.mat','simRes1_16a');
load('tempRes_16a2.mat','simRes1_16a2')
% make data arrays first
nPoints = 8;
% 
pE16a_w_simHXK_18 = ones(1,nPoints);
pE16a_w_refY3M1_18 = [1E-3 1E-2 1E-1 1E0 3E0 1E1 1E2 1E3];
% 
pE16a_error_simHXK_18 = zeros(1,nPoints);
pE16a_error_pars_18 = zeros(1,nPoints);
% 
temp_time = simRes1_16a{1}.T_FF01;
temp_hxk = simRes1_16a{1}.V_FF01(:,2);
temp_timepoints = [6:5:26, 76:50:376];
temp_HXKexp = interp1(temp_time, temp_hxk, temp_timepoints, 'pchip');
for i = 1:nPoints
    % 
    temp_HXKsim = interp1(simRes1_16a{i+1}.T_FF01, simRes1_16a{i+1}.V_FF01(:,2), temp_timepoints, 'pchip');
    pE16a_error_simHXK_18(i) = sum(abs( ( temp_HXKexp - temp_HXKsim ) ./ temp_HXKexp ));
    pE16a_error_pars_18(i) = sum(abs(xAll2(i+1,parsGLK) - x_Y3M1(parsGLK)));
end

% make data arrays second
nPoints = 16;
%
pE16a_w_simHXK_16a2 = ones(1,nPoints);
pE16a_w_refY3M1_16a2 = [3E-3,2E-2,...
                        0.1000,0.1292,0.1668,0.2154,0.2783,...
                        0.3594,0.4642,0.5995,0.7743,1.0000,...
                        3E0,3E1,3E2,3E3];
% 
pE16a2_error_simHXK = zeros(1,nPoints);
pE16a2_error_pars = zeros(1,nPoints);
% 
temp_time = simRes1_16a2{1}.T_FF01;
temp_hxk = simRes1_16a2{1}.V_FF01(:,2);
temp_timepoints = [6:5:26, 76:50:376];
temp_HXKexp = interp1(temp_time, temp_hxk, temp_timepoints, 'pchip');
for i = 1:nPoints
    %
    temp_HXKsim = interp1(simRes1_16a2{i+1}.T_FF01, simRes1_16a2{i+1}.V_FF01(:,2), temp_timepoints, 'pchip');
    pE16a2_error_simHXK(i) = sum(abs( ( temp_HXKexp - temp_HXKsim ) ./ temp_HXKexp ));
    pE16a2_error_pars(i) = sum(abs(xAll1(i+1,parsGLK) - x_Y3M1(parsGLK)));
end


%% 
% % first
% pE16a_w_simHXK_18
% pE16a_w_refY3M1_18
% pE16a_error_simHXK_18
% pE16a_error_pars_18

% % second
% pE16a_w_simHXK_16a2
% pE16a_w_refY3M1_16a2
% pE16a2_error_simHXK
% pE16a2_error_pars

% full array
full_w_simHXK = [pE16a_w_simHXK_18, pE16a_w_simHXK_16a2];
full_w_refY3M1 = [pE16a_w_refY3M1_18, pE16a_w_refY3M1_16a2];
full_error_simHXK = [pE16a_error_simHXK_18, pE16a2_error_simHXK];
full_error_pars = [pE16a_error_pars_18, pE16a2_error_pars];

% 
[temp_B,temp_I] = sort(full_w_refY3M1);
% %% % 
% figure(10001)
% semilogy(temp_B,'o-')

% 
% full_w_simHXK = [pE16a_w_simHXK, pE16a_w_simHXK_16a2]
full_w_refY3M1 = temp_B;
full_error_simHXK = full_error_simHXK(temp_I);
full_error_pars = full_error_pars(temp_I);
% % full_error_simHXK = sortrows(full_error_simHXK,full_w_refY3M1);
% % full_error_pars = sortrows(full_error_pars,full_w_refY3M1);
 
% % %%
% figure(10001)
% subplot(2,2,1), semilogy(full_w_simHXK,'o-'),title('wHXK')
% subplot(2,2,2), semilogy(full_w_refY3M1,'o-'),title('wrefY3M1')
% subplot(2,2,3), semilogx(full_w_refY3M1,full_error_simHXK,'o-'),title('simHXK')
% subplot(2,2,4), semilogx(full_w_refY3M1,full_error_pars,'o-'),title('pars')
% % subplot(2,2,3), plot(full_error_simHXK,'o-'),title('simHXK')
% % subplot(2,2,4), plot(full_error_pars,'o-'),title('pars')


%% plotting
% clear axis
% clf(1001)
figure(1001)

% semilogx
sp_h = subplot(2,2,2);
yyaxis left
% 
nVals = 100;
tempRange = linspace(-3 , 3, nVals);
tempRange_plot = logspace(-3 , 3, nVals);
% % % % temp_y = pE16a_error_simHXK_916;
temp_y = pE16a_error_simHXK_18;
temp_y(5) = [];
temp_lin = linspace(-3,3,7);
% temp_y = full_error_simHXK;
% 
% temp_lin = linspace(-3,3,7);
% temp_lin = linspace(-3,3,length(temp_y));
% 
temp_spline = spline(temp_lin, temp_y, tempRange);
    temp_spline = max(temp_spline,0.05);
    % changed spline
    [temp_lin2,ia,~] = unique(log10(full_w_refY3M1));
    temp_y2 = full_error_simHXK(ia);
    temp_spline = spline(temp_lin2, temp_y2, tempRange);
temp = semilogx(tempRange_plot, temp_spline, 'k-',... 
    'LineWidth',1);
hold on
% % 
% figure, subplot(121), plot(temp_lin), subplot(122), plot(temp_y)
% figure, subplot(121), plot(temp_lin2), subplot(122), plot(temp_y2)

% %%
% Splines right y-axis
nVals = 100;
tempRange = linspace(-3 , 3, nVals);
tempRange_plot = logspace(-3 , 3, nVals);
% % % % temp_y = pE16a_error_pars_916;
temp_y = pE16a_error_pars_18;
temp_y(5) = [];
% 
temp_lin = linspace(-3,3,7);
% 
temp_spline2 = spline(temp_lin, temp_y, tempRange);
    temp_spline2 = max(temp_spline2,0.05);
    % changed spline
    [temp_lin2,ia,~] = unique(log10(full_w_refY3M1));
    temp_y2 = full_error_pars(ia);
        temp_y2(5) = [];
        temp_lin2(5) = [];
    temp_spline2 = spline(temp_lin2, temp_y2, tempRange);
    
% 
% % % % % % % % sp_h1 = semilogx(pE16a_w_simHXK_916, pE16a_error_simHXK_916, 'ko',...
% % % % % % % %     'MarkerFaceColor','black','MarkerSize',5);
% % % % sp_h1 = semilogx(pE16a_w_refY3M1_18, pE16a_error_simHXK_18, 'ko',...
% % % %     'MarkerFaceColor','black','MarkerSize',5);
sp_h1 = semilogx(full_w_refY3M1, full_error_simHXK, 'ko',...
    'MarkerFaceColor','black','MarkerSize',5);
clear xlim
ylabel('model error ( - )')
set(gca,'YColor', 'k');
ylim([0 10])

hold on
yyaxis right
% 
temp = semilogx(tempRange_plot, temp_spline2, 'k--',...
    'LineWidth',1);
% 
% % % % % % % % % % % % % sp_h2 = semilogx(pE16a_w_simHXK_916, pE16a_error_pars_916, 'ko',...
% % % % % % % % % % % %     'MarkerSize',5);
% % % % % % % % % sp_h2 = semilogx(pE16a_w_refY3M1_18, pE16a_error_pars_18, 'ko',...
% % % % % % % % %     'MarkerSize',5);
% % % % sp_h2 = semilogx(full_w_refY3M1, full_error_pars, 'ko',...
% % % %     'MarkerSize',5);
sp_h2 = semilogx(full_w_refY3M1([1:4,7:end]), full_error_pars([1:4,7:end]), 'ko',...
    'MarkerSize',5);
ylabel('parameter deviation ( - - )')
set(gca,'YColor', 'k');%,'color','k');
xlim([0.001 1000])
ylim([0 10])

% 
tempText = {'Weighting factor'; 'parameters -> simulation'};
temp2 = text(1,-1.5,tempText,'color','black',...
    'FontSize',12,'HorizontalAlignment','center');
sp_h.Position = [0.5989    0.6130-0.025    0.2845    0.2900];


%% subplot(2,2,4) -> differente in parameter values
% clf(1001)
figure(1001)
% 
% delete(sp_h4)
sp_h4 = subplot(2,2,4);

% xFinal
% % xFinal = xAll2(15,:);
% % xFinal(30) = xAll2(14,30); 
xFinal = xAll1(11,:);

% layout
maxWidth = 5;
maxRange = 1;
nSel = 3;
colRange = cool(nSel);
xSel = xAll2([1 15 10],:);
% xSel = xAll2([1 2 5],:); % [1 2 5]
    % % <== IT IS HERE WHERE CHANGES MIGHT BE NEEDED.
    % % 
    % % 
    % % 
xSel(2,:) = xFinal;
% xSel_ref
xSel_ref = xSel;
xSel_ref(1,:) = xSel_ref(1,:) - xSel_ref(3,:); % 'normalizing' initial
xSel_ref(2,:) = xSel_ref(2,:) - xSel_ref(3,:); % 'normalizing' estimated
xSel_ref(3,:) = xSel_ref(3,:) - xSel_ref(3,:); % making it zeros

%
xlen = 1:length(parsGLK);
tempLabel = cell(length(parsGLK),1);
% clf(103)

% intial baseline -> in vitro values
area([1 7],[-1 -1],'LineStyle','none','faceColor',[0.95 0.95 0.95])
hold on
area([1 7],[1 1],'LineStyle','none','faceColor',[0.95 0.95 0.95])
line([1 7],[0 0],'Color','black','Linewidth',1)
colRange = [0.5 0.5 0.5; 0 0 1; 0 0 0];

% lines difference
for i = xlen
    %     wid = maxWidth * (abs(xSel(2,parsGLK(i)) - xSel(3,parsGLK(i)))/maxRange);
%     if i == 2
%         if xSel_ref(2,parsGLK(i)) > 0
%             line([xlen(i) xlen(i)],[0 xSel_ref(2,parsGLK(i))],...
%                 'Color','black','LineWidth',1)
%         elseif xSel_ref(2,parsGLK(i)) < 0
%             line([xlen(i) xlen(i)],[xSel_ref(2,parsGLK(i)) 0],...
%                 'Color','black','LineWidth',1)
%         end
        line([xlen(i) xlen(i)],[min([xSel_ref(2,parsGLK(i)) xSel_ref(3,parsGLK(i))]) max([xSel_ref(2,parsGLK(i)) xSel_ref(3,parsGLK(i))])],...
            'Color','black','LineWidth',1)
%     end
    hold on

    % 
    tempLabel{i} = strtok(legenda.parameters{parsGLK(i)}(9:end),', ');
end
hold on

% % plot the lines pointing at the range
for i = 1:2
    if i == 1
        plot(xlen,abs(xSel_ref(i,parsGLK)),'color','k','LineStyle','-')
        plot(xlen,-abs(xSel_ref(i,parsGLK)),'color','k','LineStyle','-')
    elseif i == 2
        plot(xlen,abs(xSel_ref(i,parsGLK)),'color','k','LineStyle','--')
        plot(xlen,-abs(xSel_ref(i,parsGLK)),'color','k','LineStyle','--')
    end
end
ylabel('reference parameter set deviation (log scale)')

% plot the parameter values in a scatter plot
% h = zeros(1,nSel);
h = zeros(1,2);
% % % % for i = 1:nSel
for i = 1:2
%     h(i) = scatter(xlen,xSel(i,parsGLK),'filled', ...
%         'MarkerEdgeColor', colRange(i,:), ...
%         'MarkerFaceColor', colRange(i,:));
    if i == 1
        h(i) = scatter(xlen,xSel_ref(i,parsGLK),'filled', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'k');
    elseif i == 2
        h(i) = scatter(xlen,xSel_ref(i,parsGLK),'filled', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', 'w');
    end
end
xlim([1 7])
xticks(xlen)
    tempLabel{1} = 'K_{m,ADP}';
    tempLabel{2} = 'K_{m,ATP}';
    tempLabel{3} = 'K_{eq}';
    tempLabel{4} = 'K_{m,G6P}';
    tempLabel{5} = 'K_{m,GLC}';
    tempLabel{6} = 'K_{i,T6P}';
    tempLabel{7} = 'V_{max}';
xticklabels(tempLabel)
xtickangle(45)
clear ylim
ylim([-3.5 3.5])
leg_h = legend(h,'non-regularized','regularized',...
    'Location','SouthOutside','Orientation','horizontal'); % ,'Box','off'
leg_h.Box = 'off';
box on
% 
sp_h4.Position = [0.5989    0.1795+0.025    0.2845    0.2900];
leg_h.Position = [0.5501+0.025 0.1392 0.3330 0.0215];


%% finalizing the figure
% make it white
% % % % set(1001,'color','w')
% % % % % box
% % % % textbox_handle = annotation('textbox');
% % % % textbox_handle.Position = [0.525 0.1 0.45 0.825];
% % % % textbox_handle.LineWidth = 2;
% text regularization
text_h = text(4,13.4,'Regularization of GLK kinetics parameters');
text_h.HorizontalAlignment = 'center';
text_h.FontSize = 12;
% arrow selected case
annotation('arrow',[0.05+0.42 0.115+0.4],[0.193 0.193],'color',cpal_broad.blue)
% % % % % labels A, B and C
% % % % text(-3, 13, 'A', 'FontSize', 20)
% % % % text(8, 13, 'B', 'FontSize', 20)
% % % % text(8, 3, 'C', 'FontSize', 20)
% arrow pointing at selected lambda
% annotation('arrow',[0.694 0.694],[0.65 0.6],'color','k')
annotation('arrow',[0.73 0.73],[0.70 0.65],'color','k')
% figh.Position = [962 42 958 954];


% % % % % %% 
% % % % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F5ABC');
% % % % % 
% % % % % %% 
% % % % % set(figh,'Units','Inches');
% % % % % pos = get(figh,'Position');
% % % % % set(figh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % % print(figh,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F5ABC','-dpdf','-r0')



%% Supplementary simulations
selResCell = simRes1_16a2([2,11]);
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;
% %% edits in the figures
close(3:4)
fh_1 = figure(1);
fh_2 = figure(2);
% %% 
fh_1.Position = [1 0.2370 0.8000 0.6917];
fh_2.Position = [1 0.2370 0.8000 0.6917];
% %% 
close all
selResCell = simRes1_16a2([2,11]);
referencePlotSimulations_manuscript


% % % % % %% 
% % % % % fh_2001 = figure(2001);
% % % % % fh_2002 = figure(2002);
% % % % % % 
% % % % % set(fh_2001,'Units','Inches');
% % % % % pos = get(fh_2001,'Position');
% % % % % set(fh_2001,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % % print(fh_2001,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS3pre','-dpdf','-r0')
% % % % % print(fh_2001,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS3pre','-dtiff','-r0')
% % % % % % 
% % % % % set(fh_2002,'Units','Inches');
% % % % % pos = get(fh_2002,'Position');
% % % % % set(fh_2002,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % % print(fh_2002,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS4pre','-dpdf','-r0')
% % % % % print(fh_2001,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS4pre','-dtiff','-r0')

