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

% load('xarray_20210129');
% load('xarray_1b'); % the one used previously in the estimations (and TREic is still depleted) 
% load('xarray_2_fix'); % good fitting of glycolysis but not of trehalose cycle (and TREic is still depleted) 
% load('xarray_36'); % there's like a balance between the fitting of glycolysis and treh cycle (and TREic is still depleted) 
% load('xarray_47'); % there's a worse balance between the fitting of glycolysis and treh cycle, but TREic is not depleted.
% load('x_Comb13_flxsglyc');
% load('x_Comb1_flx');
% load('x_Comb9_flx');
% load('x_Comb16_flx');
% load('x_Comb15_flxs_conc');
% load('x_Comb24_flxs_glyc');
% load('x_Comb21_flx');
% load('x_Comb23_flx');

% % % % load('x_Comb32_flx');

% load('xarray_X');
% load('xarray_Y');
% load('xarray_Z');

% get the paramter set from parEstFinal

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

% load('PEst_final_weights_combs_21_warray.mat')
% load('pEsts20210407_final_noNTH1_activeTvacuole\PEst_final_noNTH1_activeTvacuole_weights_combs_13_warray')
% selPars = parCombs;
% xarray(selPars) = xres;

% % % % %% updated parameter sets
% % % % load('x_Comb32_flx'); x = xarray;
% % % % load('pset_20210325.mat')
% % % % load('parSet_99b.mat')
% % % % %
% % % % parsGlyco = [28:38];
% % % % parsTre = [83:86,119:128];
% % % % parsAXP = [109:111,129:131];
% % % % parsAdd = 144:158;
% % % % %
% % % % x1 = x;
% % % % x2 = x99b;
% % % %     x2(parsGlyco) = x(parsGlyco);
% % % %     x2(parsTre) = x(parsTre);
% % % % %     x2(parsAXP) = x(parsAXP);
% % % %     x2(parsAdd) = x(parsAdd);
% % % % x3 = x2;
% % % %     x3([109 129]) = x([109 129]);
% % % % clear x, x = x3;


% folderName = 'pEsts_final_weights_combs_concs_allRuns\';
% loadName1 = 'PEst_final_weights_combs_3';
% loadName2 = '_warray.mat';
% loadName = [folderName, loadName1, loadName2];
% load(loadName)
% selPars = parCombs;
% xres = real(xres);
% xarray(selPars) = xres;

% folderName = 'pEst_final_withoutNTH1\';
% % loadName1 = 'PEst_final_weights_combs_21';
% loadName1 = 'PEst_final_weights_combs_24';
% loadName2 = '_warray.mat';
% loadName = [folderName, loadName1, loadName2];
% load(loadName)
% selPars = parCombs;
% xres = real(xres);
% xarray(selPars) = xres;



% xarray(158) = 0.0001;

% 
% loadName1 = 'PEst_final2_weights_combs_1';
% loadName2 = '_warray.mat';
% loadName = [loadName1, loadName2];
% load(loadName)
% selPars = parCombs;
% 
% xarray(selPars) = xres;

% load('x_WComb20');

% load('x_WComb35');

% xarray(158) = 2;

% xarray(154) = 1/xarray(154);
% xarray(157) = 0;
% xarray(158) = 0;
% x_original = xarray;

% xarray(157) = 0;
% xarray(158) = 0;
% % xarray(52) = 0.25;
% xarray(154) = -3.0083; %-6.0083

% from estimations:
% parsGLK = [28, 29, 30, 31, 32, 33, 34]; % k_{ADP}, k_{ATP},  k_{eq}, k_{G6P}, k_{GLC}, k_{T6P},  k_{cat},
% xres32_GLK = [ 0.8172   -0.5854    0.6283    0.2000 -0.9960   -1.1707    0.9924 ];
% xarray(parsGLK) = xres32_GLK;


% xres32 = [ 0.8172   -0.5854    0.6283    0.2000 -0.9960   -1.1707    0.9924   -0.5970 -0.5400   -0.6005    0.4253    0.4815 -0.0404   -0.6644    0.0007    0.1123  -0.7766];
% parsGLK = [28, 29, 30, 31, 32, 33, 34]; % k_{ADP}, k_{ATP},  k_{eq}, k_{G6P}, k_{GLC}, k_{T6P},  k_{cat},
% parsUGP = [144:148];
% parsTPS1 = [124:128];
% parCombs{32} = [parsGLK, parsUGP, parsTPS1];
% xarray(parCombs{32}) = xres32;

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
NumberCycles = 60;
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
% % load('pset_pE5_x1res.mat','x_pE5_x1_end'); x1 = x_pE5_x1_end;
% % load('pset_pE5_x2res.mat','x_pE5_x2_end'); x2 = x_pE5_x2_end;
% load('pset_pE5_x1res_temp.mat','x_pE5_x1_end'); x1 = x_pE5_x1_end;
% load('pset_pE5_x2res_temp.mat','x_pE5_x2_end'); x2 = x_pE5_x2_end;
load('pset_pE6_x1res.mat','x_pE6_x1_end'); x1 = x_pE6_x1_end;
load('pset_pE6_x2res.mat','x_pE6_x2_end'); x2 = x_pE6_x2_end;


%% adding a more detailed glycogen metabolism
% get experimental data,

% IC(40) = 100; % initial conncetration
setup.glycSynthDeg_separate = 1;

% %% sims
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x1,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% [T_FF01_2,Y_FF01_2,V_FF01_2] = simulate_FF(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{1}.T_FF01 = T_FF01_1;
% selResCell{1}.Y_FF01 = Y_FF01_1;
% selResCell{1}.V_FF01 = V_FF01_1;
% selResCell{2}.T_FF01 = T_FF01_2;
% selResCell{2}.Y_FF01 = Y_FF01_2;
% selResCell{2}.V_FF01 = V_FF01_2;
% 
% referencePlotSimulations
% plotMode = 0;

%%
% save('pE7_startSims.mat','selResCell')
% load('pE7_startSims.mat','selResCell')


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

    % 149_Ktre_ATHvac->163_Ktre_ATHec
    % 150_kcat_ATHvac->164_kcat_ATHec
setup.ATH_separate_EC_VAC = 1;


% %% sims
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
% disp('stopped here')


% %% Testing implementation of glycogen, as a sink
% % how the profile will look like
% expGlyc_time = dataset.FF01.fluxes_times';
% expGlyc_synth = dataset.FF01.fluxes{42}';
% expGlyc_deg = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}';
% 
% interp_expGlyc_time = 0:1:400;
% interp_expGlyc_synth = interp1(expGlyc_time,expGlyc_synth,interp_expGlyc_time,'pchip');
% interp_expGlyc_deg = interp1(expGlyc_time,expGlyc_deg,interp_expGlyc_time,'pchip');
% % interp_expGlyc_deg = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}';
% 
% figure,
% %
% subplot(1,2,1) % v_glyc_synth
% plot(interp_expGlyc_time,interp_expGlyc_synth,'k.-')
% hold on
% plot(expGlyc_time,expGlyc_synth,'bo')
% %
% subplot(1,2,2) % v_glyc_deg
% plot(interp_expGlyc_time,interp_expGlyc_deg,'k.-')
% hold on
% plot(expGlyc_time,expGlyc_deg,'bo')


%% setting the option
setup.glycogenReactionsSink = 1;
setup.dataset = dataset;
%
load('pset_pE10_xres.mat','x_pE10_start','x_pE10_end'); x = x_pE10_end;

% %% testing at start
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
% 
% %% testing the change in the mass balance
% setup.updated_bmf_Cx_ATH1ec = 1;
% 
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{2}.T_FF01 = T_FF01_1;
% selResCell{2}.Y_FF01 = Y_FF01_1;
% selResCell{2}.V_FF01 = V_FF01_1;
% 
% referencePlotSimulations
% plotMode = 0;


%% latest setup
setup.updated_bmf_Cx_ATH1ec = 1;
setup.updated_bmf_Cx_ATH1ec2 = 1;
x(166) = x(162);

% 
%%
% % sim1 = selResCell;
% % tic
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% sim2{1}.T_FF01 = T_FF01_1;
% sim2{1}.Y_FF01 = Y_FF01_1;
% sim2{1}.V_FF01 = V_FF01_1;
% 
% selResCell = sim2;
% referencePlotSimulations
% plotMode = 0;
% % toc


% % % % %% PSA on clamping
% % % % % changes in UDP_Glc incoming enzyme (UGP, TPS1)
% % % %     % UGP: kcat (x148), all (x144:148).
% % % %     % TPS1: kcat (x126), all (x124:128).
% % % % % changes in TREec-related enzymes (AGT1, ATH1)
% % % %     % AGT1: kcat (x152), all (x[151:154,158]).
% % % %     % ATH1vac: kcat (x150), all (x[149:150]).
% % % %     % ATH1ec: kcat (x164).
% % % %     
% % % % tempVal = [1 2 5 10 20,...
% % % %             30 50 100 1000 10000];
% % % % nPSAloop = 10;
% % % % xPSA = zeros(nPSAloop,166);
% % % % for o = 1:nPSAloop
% % % %     xPSA(o,:) = [x, tempVal(o)];
% % % % end
% % % %     
% % % % % xPSA = [x; x; x; x];
% % % % % xPSA(:,166) = [2 5 10 50]';
% % % % 
% % % % % simulations in parfor loop
% % % % [legenda] = legendaFull; %legenda for the names needed
% % % % metNames = legenda.metabolites;
% % % % reactNames = legenda.fluxes;
% % % % plotMode = 0;
% % % % 
% % % % % cluster = parcluster('local');
% % % % % pool = parpool(cluster,10);
% % % % % % pool = parpool(cluster,4);
% % % % % % 148
% % % % % parfor o = 1:nPSAloop
% % % % %     disp(o)
% % % % %     % simulate
% % % % %     x_sel = xPSA(o,:);
% % % % %     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % %     selResCell = cell(1,1);
% % % % %     selResCell{1}.T_FF01 = T_FF01_1;
% % % % %     selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % %     selResCell{1}.V_FF01 = V_FF01_1;
% % % % %     % save
% % % % %     saveName = sprintf('pE11d_PSA2_n%d.mat',o);
% % % % %     parsave_Y3M2_cluster_2(saveName, selResCell);
% % % % % end
% % % % % delete(pool)
% % % % % quit
% % % % 
% % % % 
% % % % %% plot PSA
% % % % % 
% % % % % xPSA
% % % % simRes1 = cell(1,nPSAloop);
% % % % for i = 1:nPSAloop
% % % %     loadName = sprintf('pE11d_PSA2_n%d.mat',i);
% % % %     load(loadName);
% % % %     simRes1(i) = selResCell;
% % % % end
% % % % % simRes1(1) = [];
% % % % % simRes1_xPSA = simRes1(1:5);
% % % % % simRes1_x152 = simRes1(6:10);
% % % % % simRes1_x164 = simRes1(11:15);
% % % % 
% % % % 
% % % % %% plot data
% % % % [legenda] = legendaFull; %legenda for the names needed
% % % % metNames = legenda.metabolites;
% % % % reactNames = legenda.fluxes;
% % % % 
% % % % choosedataset
% % % % 
% % % % % plotMode = 1; % single simulation
% % % % plotMode = 2; % multiple simulations
% % % % % plotMode = 10; % all
% % % % 
% % % % % simRes1_PSA
% % % % selResCell = simRes1;
% % % % referencePlotSimulations
% % % % 
% % % % 
% %% PSA on clamping
% % changes in UDP_Glc incoming enzyme (UGP, TPS1)
%     % UGP: kcat (x148), all (x144:148).
%     % TPS1: kcat (x126), all (x124:128).
% % changes in TREec-related enzymes (AGT1, ATH1)
%     % AGT1: kcat (x152), all (x[151:154,158]).
%     % ATH1vac: kcat (x150), all (x[149:150]).
%     % ATH1ec: kcat (x164).
%     
nPSAloop = 16;
xPSA = zeros(nPSAloop,166);
for o = 1:nPSAloop
    xPSA(o,:) = x;
end
%
xPSA(2:6,164) = xPSA(2:6,164) + [1 2 3 4 5]'; % p.ATH1_kcat_ec
%
xPSA(7:11,152) = xPSA(7:11,152) + [1 2 3 4 5]'; % p.AGT1_kcat
%
xPSA(12:16,164) = xPSA(12:16,164) + [1 2 3 4 5]'; % p.ATH1_kcat_ec
xPSA(12:16,152) = xPSA(12:16,152) + [1 2 3 4 5]'; % p.AGT1_kcat

xPSA_pre = xPSA;

% xPSA = [x; x; x; x];
% xPSA(:,166) = [2 5 10 50]';

% simulations in parfor loop
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
plotMode = 0;

% % % % % cluster = parcluster('local');
% % % % % pool = parpool(cluster,16);
% % % % % pool = parpool(cluster,4);
% % % % % 148
% % % % % parfor o = 1:nPSAloop
% % % % for o = 7:11
% % % %     disp(o)
% % % %     % simulate
% % % %     x_sel = xPSA(o,:);
% % % %     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     selResCell = cell(1,1);
% % % %     selResCell{1}.T_FF01 = T_FF01_1;
% % % %     selResCell{1}.Y_FF01 = Y_FF01_1;
% % % %     selResCell{1}.V_FF01 = V_FF01_1;
% % % %     % save
% % % % %     saveName = sprintf('pE11d_PSA3_n%d.mat',o);
% % % %     saveName = sprintf('pE11d_PSA3b_num1_n%d.mat',o);
% % % %     parsave_Y3M2_cluster_2(saveName, selResCell);
% % % % end
% % % % % delete(pool)
% % % % % quit
% % % % 
% % % % 
% % % % %% plot PSA
% % % % % 
% % % % % xPSA
% % % % simRes1 = cell(1,nPSAloop);
% % % % % for i = 1:nPSAloop
% % % % for i = 7:11
% % % % %     loadName = sprintf('pE11d_PSA3_n%d.mat',i);
% % % %     loadName = sprintf('pE11d_PSA3b_num1_n%d.mat',i);
% % % % %     loadName = sprintf('pE11d_PSA3b_n%d.mat',i);
% % % %     load(loadName);
% % % %     simRes1(i) = selResCell;
% % % % end
% % % % % % simRes1(1) = [];
% % % % % simRes1_x164 = simRes1([1,2:6]);
% % % % % simRes1_x152 = simRes1([1,7:11]);
% % % % % simRes1_xPSA = simRes1([1,12:16]);
% % % % 
% % % % 
% % % % %% plot data
% % % % [legenda] = legendaFull; %legenda for the names needed
% % % % metNames = legenda.metabolites;
% % % % reactNames = legenda.fluxes;
% % % % 
% % % % choosedataset
% % % % 
% % % % % plotMode = 1; % single simulation
% % % % plotMode = 2; % multiple simulations
% % % % % plotMode = 10; % all
% % % % 
% % % % % % simRes1_PSA
% % % % % selResCell = simRes1_x164;
% % % % % referencePlotSimulations
% % % % % 
% % % % % % simRes1_PSA
% % % % % selResCell = simRes1_x152;
% % % % % referencePlotSimulations
% % % % % 
% % % % % % simRes1_PSA
% % % % % selResCell = simRes1_xPSA;
% % % % % referencePlotSimulations
% % % % 
% % % % selResCell = simRes1(7:11);
% % % % referencePlotSimulations


% % % % % selResCell([1:6,12:16]) = [];
% % % % % referencePlotSimulations
% 
% % disp('something');
% 
% %% PSA on clamping
% % changes in UDP_Glc incoming enzyme (UGP, TPS1)
%     % UGP: kcat (x148), all (x144:148).
%     % TPS1: kcat (x126), all (x124:128).
% % changes in TREec-related enzymes (AGT1, ATH1)
%     % AGT1: kcat (x152), all (x[151:154,158]).
%     % ATH1vac: kcat (x150), all (x[149:150]).
%     % ATH1ec: kcat (x164).
% parsATH1ec = 163:165;
% parsAGT1 = [151:154,158,166];
% lenPars = length([parsAGT1, parsATH1ec]);
% 
% nPSAloop = 16*10;
% xPSA = zeros(nPSAloop,166);
% for o = 1:nPSAloop
%     xPSA(o,:) = x;
% end
% rng(1), xPSA(:,[parsAGT1, parsATH1ec]) = ...
%     xPSA(:,[parsAGT1, parsATH1ec]) + 10*rand(nPSAloop,lenPars);
% 
% 
% % xPSA = [x; x; x; x];
% % xPSA(:,166) = [2 5 10 50]';
% 
% % simulations in parfor loop
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% plotMode = 0;
% 
% % cluster = parcluster('local');
% % pool = parpool(cluster,16);
% % % pool = parpool(cluster,4);
% % % 148
% % parfor o = 1:nPSAloop
% %     disp(o)
% %     % simulate
% %     x_sel = xPSA(o,:);
% %     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% %     selResCell = cell(1,1);
% %     selResCell{1}.T_FF01 = T_FF01_1;
% %     selResCell{1}.Y_FF01 = Y_FF01_1;
% %     selResCell{1}.V_FF01 = V_FF01_1;
% %     % save
% %     saveName = sprintf('pE11d_PSA4_n%d.mat',o);
% %     parsave_Y3M2_cluster_2(saveName, selResCell);
% % end
% % delete(pool)
% % quit
% 
% 
% %% plot PSA
% % 
% % xPSA
% simRes1 = cell(1,nPSAloop);
% for i = 1:nPSAloop
%     loadName = sprintf('pE11d_PSA4_n%d.mat',i);
%     load(loadName);
%     simRes1(i) = selResCell;
% end
% % % simRes1(1) = [];
% % simRes1_x164 = simRes1([1,2:6]);
% % simRes1_x152 = simRes1([1,7:11]);
% % simRes1_xPSA = simRes1([1,12:16]);
% 
% 
% %% plot data
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% 
% choosedataset
% 
% % plotMode = 1; % single simulation
% plotMode = 2; % multiple simulations
% % plotMode = 10; % all
% 
% % simRes1_PSA
% selResCell = simRes1;
% referencePlotSimulations
% 
% 
% % %% Get the clamped setup again
% % setup.clampd_ATHec_AGT1 = 1;
% % setup.clamp_factor_ATHec = 10;
% % 
% % [legenda] = legendaFull; %legenda for the names needed
% % metNames = legenda.metabolites;
% % reactNames = legenda.fluxes;
% % % plotflag = 2; % variables by iteration
% % choosedataset
% % setup.experiment = 1;
% % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % plotMode = 2; % multiple simulations
% % % clear selResCell
% % selResCell{1}.T_FF01 = T_FF01_1;
% % selResCell{1}.Y_FF01 = Y_FF01_1;
% % selResCell{1}.V_FF01 = V_FF01_1;
% % % selResCell(1) = []
% % referencePlotSimulations
% % plotMode = 0;
% % 
% % setup.clampd_ATHec_AGT1 = 0;
% 
% 
% %% Another way
% % sim1 = selResCell;
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% sim2{1}.T_FF01 = T_FF01_1;
% sim2{1}.Y_FF01 = Y_FF01_1;
% sim2{1}.V_FF01 = V_FF01_1;
% 
% selResCell = sim2;
% referencePlotSimulations
% plotMode = 0;
% 
% 
% %% checking equation
% 
% % ATH1ec
% f.ATH1ec = 0.00196 * 0.9;
% p.ATH1_kcat_ec = 10 .^ x(164) .* 99 / 60 / 2e-3 / 0.00196;
% p.ATH1_Ktre_ec = 10 .^ x(163) .* 4.7 ; 
% p.ATH1_Kt6p_ec = 10 .^ x(165) .* 0.1; %(UNIT!)
% TREec = sim2{1}.Y_FF01(:,37);
% T6P = sim2{1}.Y_FF01(:,26);
% v_ATH1ec_calc = (p.ATH1_kcat_ec .* f.ATH1ec).* ( (TREec./p.ATH1_Ktre_ec)./ ( 1+(TREec./p.ATH1_Ktre_ec)+(T6P./p.ATH1_Kt6p_ec) ) );
% v_ATH1ec_sim = sim2{1}.V_FF01(:,45);
% v_ATH1ec_clamped = interp1(setup.dataset.FF01.fluxes_times', 10*0.9*setup.dataset.FF01.fluxes{34}', sim2{1}.T_FF01, 'pchip');
% 
% % AGT1
% f.AGT1 = 0.000067;
% p.AGT1_kcat = 10.^x(152).* 0.140 / 60 / 2e-3 / 0.000067;
% p.AGT1_Ktre = 10.^x(153).* 4;
% p.AGT1_Ktre_ec = 10.^x(151).* 4;
% p.AGT1_Kt6p = 10 .^ x(166) .* 0.1;
% p.AGT1_Keq = 10.^x(154).* 1000;
% TRE = sim2{1}.Y_FF01(:,25);
% TREec = sim2{1}.Y_FF01(:,37);
% T6P = sim2{1}.Y_FF01(:,26);
% v_AGT1_calc = (p.AGT1_kcat .* f.AGT1).* (1./ p.AGT1_Ktre) .* ( TRE - TREec./p.AGT1_Keq )./ ( 1+TRE./p.AGT1_Ktre + TREec./p.AGT1_Ktre_ec + (T6P./p.AGT1_Kt6p) ) ;
% v_AGT1_sim = sim2{1}.V_FF01(:,47);
% v_AGT1_clamped = interp1(setup.dataset.FF01.fluxes_times', dataset.FF01.fluxes{32}' - dataset.FF01.fluxes{33}', sim2{1}.T_FF01, 'pchip');
% 
% figure,
% subplot(2,3,1), plot(sim2{1}.T_FF01,v_ATH1ec_calc,'o-'), title('v_ATH1ec_calc')
% subplot(2,3,2), plot(sim2{1}.T_FF01,v_ATH1ec_sim,'o-'), title('v_ATH1ec_sim')
% subplot(2,3,3), plot(sim2{1}.T_FF01,v_ATH1ec_clamped,'o-'), title('v_ATH1ec_clamped')
% subplot(2,3,4), plot(sim2{1}.T_FF01,v_AGT1_calc,'o-'), title('v_AGT1_calc')
% subplot(2,3,5), plot(sim2{1}.T_FF01,v_AGT1_sim,'o-'), title('v_AGT1_sim')
% subplot(2,3,6), plot(sim2{1}.T_FF01,v_AGT1_clamped,'o-'), title('v_AGT1_clamped')


%% PSA on clamping
% changes in UDP_Glc incoming enzyme (UGP, TPS1)
    % UGP: kcat (x148), all (x144:148).
    % TPS1: kcat (x126), all (x124:128).
% changes in TREec-related enzymes (AGT1, ATH1)
    % AGT1: kcat (x152), all (x[151:154,158]).
    % ATH1vac: kcat (x150), all (x[149:150]).
    % ATH1ec: kcat (x164).
    
nPSAloop = 31;
xPSA = zeros(nPSAloop,166);
for o = 1:nPSAloop
    xPSA(o,:) = x;
end
%
xPSA(1:5,151) = xPSA(1:5,151) + [1 2 3 4 5]'; % p.ATH1_kcat_ec
    xPSA(6:10,152) = xPSA(6:10,152) + [1 2 3 4 5]'; % p.AGT1_kcat
xPSA(11:15,153) = xPSA(11:15,153) + [1 2 3 4 5]'; % p.AGT1_Ktre
xPSA(16:20,154) = xPSA(16:20,154) + [1 2 3 4 5]'; % p.AGT1_Keq
xPSA(21:25,158) = xPSA(21:25,158) + [1 2 3 4 5]'; % p.AGT1_Ki
xPSA(26:30,166) = xPSA(26:30,166) + [1 2 3 4 5]'; % p.AGT1_Ki

% p.AGT1_Ktre_ec  = 10.^x(151).* 4; % mM (Stambuk et al.1998)
% p.AGT1_kcat     = 10.^x(152).* AGT1_kcat;
% p.AGT1_Ktre     = 10.^x(153).* 4; % mM (Stambuk et al.1998)
% p.AGT1_Keq      = 10.^x(154).* 1000; 
% p.AGT1_Ki       = 10.^x(158).* 1; % negative correlation with UDPG

% simulations in parfor loop
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
plotMode = 0;


% cluster = parcluster('local');
% pool = parpool(cluster,4);
% % pool = parpool(cluster,16);
% parfor o = 1:nPSAloop
% % for o = 6:10
%     disp(o)
%     % simulate
%     x_sel = xPSA(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE11d_PSA10_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% delete(pool)
% quit


%% plot PSA
% 
% xPSA
simRes1 = cell(1,nPSAloop);
for i = 1:nPSAloop
% for i = 6:10
%     loadName = sprintf('pE11d_PSA5c_num2_n%d.mat',i);
    loadName = sprintf('pE11d_PSA10_n%d.mat',i);
    load(loadName);
    simRes1(i) = selResCell;
end
% simRes1(1) = [];
simRes1_x151 = simRes1([1:5,31]);
simRes1_x152 = simRes1([6:10,31]);
simRes1_x153 = simRes1([11:15,31]);
simRes1_x154 = simRes1([16:20,31]);
simRes1_x158 = simRes1([21:25,31]);
simRes1_x166 = simRes1([26:30,31]);
% xPSA(1:5,151) = xPSA(1:5,151) + [1 2 3 4 5]'; % p.ATH1_kcat_ec
% xPSA(6:10,152) = xPSA(6:10,152) + [1 2 3 4 5]'; % p.AGT1_kcat
% xPSA(11:15,153) = xPSA(11:15,153) + [1 2 3 4 5]'; % p.AGT1_Ktre
% xPSA(16:20,154) = xPSA(16:20,154) + [1 2 3 4 5]'; % p.AGT1_Keq
% xPSA(21:25,158) = xPSA(21:25,158) + [1 2 3 4 5]'; % p.AGT1_Ki
% xPSA(26:30,158) = xPSA(26:30,166) + [1 2 3 4 5]'; % p.AGT1_Ki


%% plot data
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;

choosedataset

% plotMode = 1; % single simulation
plotMode = 2; % multiple simulations
% plotMode = 10; % all

xPSA(1:5,151) = xPSA(1:5,151) + [1 2 3 4 5]'; % p.ATH1_kcat_ec
xPSA(6:10,152) = xPSA(6:10,152) + [1 2 3 4 5]'; % p.AGT1_kcat
xPSA(11:15,153) = xPSA(11:15,153) + [1 2 3 4 5]'; % p.AGT1_Ktre
xPSA(16:20,154) = xPSA(16:20,154) + [1 2 3 4 5]'; % p.AGT1_Keq
xPSA(21:25,158) = xPSA(21:25,158) + [1 2 3 4 5]'; % p.AGT1_Ki
xPSA(26:30,166) = xPSA(26:30,166) + [1 2 3 4 5]'; % p.AGT1_Ki

% simRes1_PSA151
selResCell = simRes1_x151;
referencePlotSimulations

% simRes1_PSA152
selResCell = simRes1_x152;
referencePlotSimulations

% simRes1_PSA153
selResCell = simRes1_x153;
referencePlotSimulations

% simRes1_PSA154
selResCell = simRes1_x154;
referencePlotSimulations

% simRes1_PSA158
selResCell = simRes1_x158;
referencePlotSimulations

% simRes1_PSA166
selResCell = simRes1_x166;
referencePlotSimulations

