% % PAREST1_REPLICATE_INITIAL.m
% In this code the tre cycle and glt, hxk are optimized as it happen before
% to obtain the parameter set x32.

% structure is:
    % early sampling
    % parameter estimation with paraemter and cost function combinations
    

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
load('pset_pE2res.mat');
x1 = x_pE2_start;
x2 = x_pE2_end;
%% glycerol_SynthDeg active
setup.glycSynthDeg = 1; % this adds a reaction to the mass balance (ODEs), reaction calculation in rate equations and the starting parameter value
x1(159) = -4; % -10;
x2(159) = -4; %-10;
%% extension glycogen metabolism
% load parameters
load('pset_pE6_x1res.mat','x_pE6_x1_end'); x1 = x_pE6_x1_end;
load('pset_pE6_x2res.mat','x_pE6_x2_end'); x2 = x_pE6_x2_end;
%% adding a more detailed glycogen metabolism
setup.glycSynthDeg_separate = 1;
%% recall last case
x = x1;
x(162) = 10;
setup.ATHinhibitionT6P = 1;
load('pset_pE7_xres.mat','x_pE7_start','x_pE7_end');
clear x x1 x2
x = x_pE7_end;
%% separating ath1_ec and ath1_vac
x(163:164) = x(149:150); % km_tre, k_cat
x(165) = x(162); % ki_t6p
setup.ATH_separate_EC_VAC = 1;
%% setting the option
setup.glycogenReactionsSink = 1;
setup.dataset = dataset;
%
load('pset_pE10_xres.mat','x_pE10_start','x_pE10_end'); x = x_pE10_end;
%% latest setup
% % % % setup.updated_bmf_Cx_ATH1ec = 1;
setup.TREec_brothOut_OFF = 1;
setup.updated_bmf_Cx_ATH1ec = 1;


%% added for enrichment simulations
% % % % load('datasetEnrich.mat');
% % % % reorganiseEnrichData;


% %% simulation enrichment
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% legendaMetabolites_addEnrichment;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{1}.T_FF01 = T_FF01_1;
% selResCell{1}.Y_FF01 = Y_FF01_1;
% selResCell{1}.V_FF01 = V_FF01_1;
% % 
% referencePlotSimulations_enrichment
% plotMode = 0;


% %% development of trehalose secretion: clamping incoming glc_ec
% % Start with change
% % setup.clamp_GLCec = 0;
% setup.clamp_GLCec = 1;
% %
% NumberCycles = 5;
% 
% % %% testing the change in the mass balance
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{1}.T_FF01 = T_FF01_1;
% selResCell{1}.Y_FF01 = Y_FF01_1;
% selResCell{1}.V_FF01 = V_FF01_1;
% 
% referencePlotSimulations
% plotMode = 0;
% 
% setup.clamp_GLCec = 0;

%% parameter estimation setup
% latest added
setup.clamp_GLCec = 1;
NumberCycles = 5;

% %% testing the change in the mass balance
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
% plotflag = 2; % variables by iteration
choosedataset


% % from here pEst setup

% % previous setup
% setup.changing_Keq_glt = 1;
% tempValues = -[0 0.1 0.5 1 2 3];
% setup.Keq_glt_inc = tempValues(i);
% setup.clamp_GLCec = 1;

setup.changing_Keq_glt = 2; % here the value number 3.
x(166) = 0; % % <== change x166 if needed here.


% %% testing the change in the mass balance
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
% plotflag = 2; % variables by iteration
choosedataset
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
% clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;

referencePlotSimulations
plotMode = 0;

setup.clamp_GLCec = 0;


%% check ex-vivo the functions proposed in CA_SuarezMendez
% exp data
S_data = interp1(dataset.FF01.metabolites.ECglucose.time(1:end-2),... 
    dataset.FF01.metabolites.ECglucose.conc(1:end-2),...
    dataset.FF01.fluxes_times,'pchip'); % mM
Vglt_data = dataset.FF01.fluxes{1}; % mM s-1
% simpleMM approach
CASM_simpleMM.qsmax = 279; % mmol Cmol-1 h-1
CASM_simpleMM.Km = 1; % mM
CASM_simpleMM.Vglt_sim = CASM_simpleMM.qsmax .* S_data ./ (CASM_simpleMM.Km + S_data); % mmol Cmol-1 h-1
% minGluc approach
CASM_minGluc.Csmin = 0.05; % mM
CASM_minGluc.Csmin_exp = S_data(1); %0.05; % mM
CASM_minGluc.K = 214.6; % L Cmol-1 h-1
CASM_minGluc.Vglt_sim = CASM_minGluc.K .* (S_data - CASM_minGluc.Csmin) ; % mmol Cmol-1 h-1
CASM_minGluc.Vglt_sim_data = CASM_minGluc.K .* (S_data - CASM_minGluc.Csmin_exp) ; % mmol Cmol-1 h-1

% check the units (to accurate yet)
ConvFactor = 0.7240/85; % from (mmol Cmol-1 h-1) to (mM s-1)
CASM_simpleMM.Vglt_sim_units = CASM_simpleMM.Vglt_sim * ConvFactor;
CASM_minGluc.Vglt_sim_units = CASM_minGluc.Vglt_sim * ConvFactor;
CASM_minGluc.Vglt_sim_data_units = CASM_minGluc.Vglt_sim_data * ConvFactor;

%%
% plotting
clf(10)
figure(10),
plot(dataset.FF01.fluxes_times, CASM_simpleMM.Vglt_sim_units,'k.-')
hold on
plot(dataset.FF01.fluxes_times, CASM_minGluc.Vglt_sim_units,'k.--')
hold on
plot(dataset.FF01.fluxes_times, CASM_minGluc.Vglt_sim_data_units,'k.:')
hold on
plot(dataset.FF01.fluxes_times, Vglt_data,'ko','MarkerFaceColor','r')
legend('simple MM.', 'min. Glucose 0.05', 'min. Glucose 0.1', 'exp.')
xlabel('time (s)')
ylabel('Vglt (mM s^{-1})')


%% testing this in the full model

% blank and constant setup
blankWeight = zeros(1,88); % 85+3 for glycerol
    warray = blankWeight;
lambdalist = 0;
setup.parEst.lambda = lambdalist(1); lam = setup.parEst.lambda;
options = optimoptions('lsqnonlin','Display','iter',...
    'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4,...
    'DiffMinChange',0.1);

% core options:
x_temp = zeros(1,2); %xArray(o,selPars);
plength = length(x_temp);

% boundaries
lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;

% parameter estimation
tic
[xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGLT_CASM_minGluc,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,NumberCycles,IC,selPars,warray);
t = toc; 


%% outcoming simulations
% minGluc approach
EstModel_minGluc.Csmin = 0.05 .* 10 .^ xres(1); % mM
EstModel_minGluc.K = 214.6 .* 10 .^ xres(2); % L Cmol-1 h-1
EstModel_minGluc.Vglt_sim = EstModel_minGluc.K .* (S_data - EstModel_minGluc.Csmin) ; % mmol Cmol-1 h-1

% check the units (to accurate yet)
ConvFactor = 0.7240/85; % from (mmol Cmol-1 h-1) to (mM s-1)
EstModel_minGluc.Vglt_sim_units = EstModel_minGluc.Vglt_sim * ConvFactor;

% plotting
% clf(11)
figure(11),
plot(dataset.FF01.fluxes_times, EstModel_minGluc.Vglt_sim_units,'k.-')
hold on
plot(dataset.FF01.fluxes_times, CASM_minGluc.Vglt_sim_units,'k.--')
hold on
plot(dataset.FF01.fluxes_times, CASM_minGluc.Vglt_sim_data_units,'k.:')
hold on
plot(dataset.FF01.fluxes_times, Vglt_data,'ko','MarkerFaceColor','r')
legend('Estimated model', 'min. Glucose 0.05', 'min. Glucose 0.1', 'exp.')
xlabel('time (s)')
ylabel('Vglt (mM s^{-1})')


%% Directly implementing the Csmin idea in out model
% latest added
setup.clamp_GLCec = 1;
% setup.clamp_GLCec = 0;
NumberCycles = 5;

% %% testing the change in the mass balance
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
% plotflag = 2; % variables by iteration
choosedataset


% % from here pEst setup

% % previous setup
% setup.changing_Keq_glt = 1;
% tempValues = -[0 0.1 0.5 1 2 3];
% setup.Keq_glt_inc = tempValues(i);
% setup.clamp_GLCec = 1;

% setup.changing_Keq_glt = 2; % here the value number 3.
setup.changing_Keq_glt = 3; % here the value number 3.
% setup.changing_Keq_glt = 0; % here the value number 3.
x(166) = 0; % % <== change x166 if needed here.


% %% testing the change in the mass balance
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
% plotflag = 2; % variables by iteration
choosedataset
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
% clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;

referencePlotSimulations
plotMode = 0;

setup.clamp_GLCec = 0;


