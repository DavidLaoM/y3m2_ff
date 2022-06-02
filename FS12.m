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
if setup.GPdataset.GP400WT == 1
    InitCond_GP_TUD;
elseif setup.GPdataset.GP1800WT == 1
    InitCond_GP_TUD_1800;
elseif setup.GPdataset.GP400M == 1
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


%% (from pEst18a) add final change (Km or pSet)
% changes in parameters to account for only ATPase change
load('pset_Y3M1.mat'); 
x_Y3M1 = x99b;
x(109) = x_Y3M1(109); % ATPase changed here directly
x(129) = x_Y3M1(129) + 1; % ATPase changed here directly
% settin the same Km value (km_i = km_o)
x2 = x; x2(35) = x2(38);

% % % % % %%
% % % % % load('parSet_105b.mat','x105b');
% % % % % x2 = x;
% % % % %     parsPFK = 43:56;
% % % % % %     parsALD = 11:15;
% % % % %     parsPYK = 71:77;
% % % % % x2(parsPFK) = x105b(parsPFK);
% % % % % % x2(parsALD) = x105b(parsALD);
% % % % % x2(parsPYK) = x105b(parsPYK);
% % % % % x2(152) = -2.4128;

%% randomize the parameter set
N = 10000; %test multiple combinations
% rng(1), rand_dist = rand(N,length(x2)) * 0.2 - 0.1;
rng(1), rand_dist = rand(N,length(x2)) * (log10(1.1)-log10(0.9)) + log10(0.9);
% 20 + (150-20) 

% figure, histogram(rand_dist)
xRobustness = zeros(N,length(x2));
for i = 1:N
    xRobustness(i,:) = x2 + rand_dist(i,:);
end


%% Simulation reference FF01
% NumberCycles = 5;
% % first simulation with the reference case
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% selResCell{1}.T_FF01 = T_FF01_1;
% selResCell{1}.Y_FF01 = Y_FF01_1;
% selResCell{1}.V_FF01 = V_FF01_1;
% % loop to simulate the robustness check
% simRobust = cell(1,N);
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% parfor i = 1:N
%     disp(i)
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(xRobustness(i,:),canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRobust{i}.T_FF01 = T_FF01_1;
%     simRobust{i}.Y_FF01 = Y_FF01_1;
%     simRobust{i}.V_FF01 = V_FF01_1;
% end
% save('simRes_mF8_robust.mat','selResCell','simRobust')
% delete(pool)
% 
% % 
% NumberCycles = 5;
% 
% % first simulation with the reference case
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% selResIntial{1}.T_FF01 = T_FF01_1;
% selResIntial{1}.Y_FF01 = Y_FF01_1;
% selResIntial{1}.V_FF01 = V_FF01_1;
% 
% % loop to simulate the robustness check
% simRobust_r1 = cell(1,2500);
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% parfor i = 1:2500
%     disp(i)
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(xRobustness(i,:),canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRobust_r1{i}.T_FF01 = T_FF01_1;
%     simRobust_r1{i}.Y_FF01 = Y_FF01_1;
%     simRobust_r1{i}.V_FF01 = V_FF01_1;
% end
% save('simRes_mF8_robust_round1.mat','selResIntial','simRobust_r1')
% delete(pool)
% 
% % loop to simulate the robustness check
% simRobust_r2 = cell(1,2500);
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% parfor i = 2501:5000
%     disp(i)
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(xRobustness(i,:),canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRobust_r2{i-2500}.T_FF01 = T_FF01_1;
%     simRobust_r2{i-2500}.Y_FF01 = Y_FF01_1;
%     simRobust_r2{i-2500}.V_FF01 = V_FF01_1;
% end
% save('simRes_mF8_robust_round2.mat','selResIntial','simRobust_r2')
% delete(pool)
% 
% % loop to simulate the robustness check
% simRobust_r3 = cell(1,2500);
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% parfor i = 5001:7500
%     disp(i)
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(xRobustness(i,:),canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRobust_r3{i-5000}.T_FF01 = T_FF01_1;
%     simRobust_r3{i-5000}.Y_FF01 = Y_FF01_1;
%     simRobust_r3{i-5000}.V_FF01 = V_FF01_1;
% end
% save('simRes_mF8_robust_round3.mat','selResIntial','simRobust_r3')
% delete(pool)
% 
% % loop to simulate the robustness check
% simRobust_r4 = cell(1,2500);
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% parfor i = 7501:10000
%     disp(i)
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(xRobustness(i,:),canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRobust_r4{i-7500}.T_FF01 = T_FF01_1;
%     simRobust_r4{i-7500}.Y_FF01 = Y_FF01_1;
%     simRobust_r4{i-7500}.V_FF01 = V_FF01_1;
% end
% save('simRes_mF8_robust_round4.mat','selResIntial','simRobust_r4')
% delete(pool)


% %% load and better split results array
% load('simRes_mF8_robust.mat','selResCell','simRobust')
% selResIntial = selResCell;


%% Load the second one
% 
load('simRes_mF8_robust_round1.mat','selResIntial','simRobust_r1')
load('simRes_mF8_robust_round2.mat','simRobust_r2')
load('simRes_mF8_robust_round3.mat','simRobust_r3')
load('simRes_mF8_robust_round4.mat','simRobust_r4')
simRobust = [simRobust_r1, simRobust_r2, simRobust_r3, simRobust_r4];
clear simRobust_r1 simRobust_r2 simRobust_r3 simRobust_r4

%%
selResCell = [selResIntial, simRobust(1:10)];
plotMode = 2;
% referencePlotSimulations_enrichment
referencePlotSimulations
% plot_robustnessSimulations
plotMode = 0;
clear selResCell


%% visualize
plot_robustnessSimulations


% % % % % %% saving
% % % % % % fig_h1
% % % % % set(fig_h,'Units','Inches');
% % % % % pos = get(fig_h,'Position');
% % % % % set(fig_h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % % print(fig_h,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS1pre','-dpdf','-r0')
% % % % % print(fig_h,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS1pre','-dtiff','-r0')
% % % % % % fig_h2
% % % % % set(fig_h2,'Units','Inches');
% % % % % pos = get(fig_h2,'Position');
% % % % % set(fig_h2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % % print(fig_h2,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS2pre','-dpdf','-r0')
% % % % % print(fig_h2,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS2pre','-dtiff','-r0')
% % % % % %%
% % % % % savefig(fig_h,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS1pre')
% % % % % savefig(fig_h2,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\FS2pre')

