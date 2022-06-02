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
load('datasetEnrich.mat');
reorganiseEnrichData;


%% simulation
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
legendaMetabolites_addEnrichment;
% plotflag = 2; % variables by iteration
choosedataset
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
% clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;
% 
referencePlotSimulations_enrichment
plotMode = 0;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Manual development starts from here: % % % % % % % % % % % % % % % % % % 
% % legenda.metabolites{65}
% % ExpData.metabolites{65}
% temp = ExpData.metabolites{65};
% 
% figure,
% subplot(1,2,1), plot(temp.time, temp.fraction, '.-'), title('fraction')
% subplot(1,2,2), plot(temp.conc_time, temp.conc, '.-'), title('conc')
% % subplot(1,3,3), plot(temp.time, temp.conc_time, '.-'), title('fraction')

setup.clamp_enrichment_GLCec = 0;
setup.clamp_enrichment_GLCec_data = ExpData;

[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
legendaMetabolites_addEnrichment;
% plotflag = 2; % variables by iteration
choosedataset
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
% clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;
% 
referencePlotSimulations_enrichment
plotMode = 0;

setup.clamp_enrichment_GLCec = 0;


%% development of trehalose secretion: simple idea for secretion from vacuole
setup.TRE_secretion_dev = 1;

[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
legendaMetabolites_addEnrichment;
% plotflag = 2; % variables by iteration
choosedataset
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
% clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;
% 
referencePlotSimulations_enrichment
plotMode = 0;
% 
setup.TRE_secretion_dev = 0;


%% development of trehalose secretion: clamping incoming glc_ec
% Start with change
% setup.clamp_GLCec = 0;
setup.clamp_GLCec = 1;
%
NumberCycles = 5;

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


%% then a PSA on km_GLCin_glk + clamping
% Start with change + PSA on km_GLCin_glk
setup.clamp_GLCec = 1;%
NumberCycles = 5;

xPSA32 = [x; x; x; x; x];
    xPSA32(2,32) = xPSA32(1,32) + 0.25;
    xPSA32(3,32) = xPSA32(1,32) + 1;
    xPSA32(4,32) = xPSA32(1,32) + 2;
    xPSA32(5,32) = xPSA32(1,32) + 3;
    
xPSA32_2 = [x; x; x; x; x];
    xPSA32_2(2,32) = xPSA32_2(1,32) - 0.25;
    xPSA32_2(3,32) = xPSA32_2(1,32) - 1;
    xPSA32_2(4,32) = xPSA32_2(1,32) - 2;
    xPSA32_2(5,32) = xPSA32_2(1,32) - 3;
    
% simulations in parfor loop
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
plotMode = 0;
%
choosedataset
setup.experiment = 1;

% cluster = parcluster('local');
% % pool = parpool(cluster,16);
% pool = parpool(cluster,4);
% % 164
% parfor o = 1:5
%     % simulate
%     x_sel = xPSA32(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
% %     saveName = sprintf('pE11h_psa_x164_n%d.mat',o);
% %     saveName = sprintf('pE11h_psa_x164_n%d_oldSetup.mat',o);
%     saveName = sprintf('pE12a_temp_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% %
% setup.clamp_GLCec = 0;
% 
% % delete(pool)
% quit


%%
% 164
for o = 1:5
    % simulate
%     x_sel = xPSA32(o,:);
    x_sel = xPSA32_2(o,:);
    [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    selResCell = cell(1,1);
    selResCell{1}.T_FF01 = T_FF01_1;
    selResCell{1}.Y_FF01 = Y_FF01_1;
    selResCell{1}.V_FF01 = V_FF01_1;
    % save
    saveName = sprintf('pE12a_temp2_n%d.mat',o);
    parsave_Y3M2_cluster_2(saveName, selResCell);
end


%%
%
simRes1 = cell(1,5);
for i = 1:5
%     loadName = sprintf('pE12a_temp_n%d.mat',i);
    loadName = sprintf('pE12a_temp2_n%d.mat',i);
    load(loadName);
    simRes1(i) = selResCell;
end
% simRes1
plotMode = 2; % multiple simulations
selResCell = simRes1;
referencePlotSimulations


%% Then we try this for glt
xBlank = [x; x; x; x];
% 35
xPSA35up = xBlank; xPSA35up(:,35) = xPSA35up(:,35) + [0 .5 1 2]';
xPSA35down = xBlank; xPSA35down(:,35) = xPSA35down(:,35) - [0 .5 1 2]';
% 36
xPSA36up = xBlank; xPSA36up(:,36) = xPSA36up(:,36) + [0 .5 1 2]';
xPSA36down = xBlank; xPSA36down(:,36) = xPSA36down(:,36) - [0 .5 1 2]';
% 38
xPSA38up = xBlank; xPSA38up(:,38) = xPSA38up(:,38) + [0 .5 1 2]';
xPSA38down = xBlank; xPSA38down(:,38) = xPSA38down(:,38) - [0 .5 1 2]';

% % parallel
% cluster = parcluster('local');
% pool = parpool(cluster,4);
% % 35up
% parfor o = 1:4
%     % simulate
%     x_sel = xPSA35up(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE12a_temp35u_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% % 35down
% parfor o = 1:4
%     % simulate
%     x_sel = xPSA35down(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE12a_temp35d_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% % 36up
% parfor o = 1:4
%     % simulate
%     x_sel = xPSA36up(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE12a_temp36u_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% % 36down
% parfor o = 1:4
%     % simulate
%     x_sel = xPSA36down(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE12a_temp36d_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% % 38up
% parfor o = 1:4
%     % simulate
%     x_sel = xPSA38up(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE12a_temp38u_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% % 38down
% parfor o = 1:4
%     % simulate
%     x_sel = xPSA38down(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE12a_temp38d_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
% end
% % selResCell


%%
%
simRes1 = cell(1,4);
for i = 1:4
%     loadName = sprintf('pE12a_temp35u_n%d.mat',i);
%     loadName = sprintf('pE12a_temp35d_n%d.mat',i);
%     loadName = sprintf('pE12a_temp36u_n%d.mat',i);
%     loadName = sprintf('pE12a_temp36d_n%d.mat',i);
%     loadName = sprintf('pE12a_temp38u_n%d.mat',i);
    loadName = sprintf('pE12a_temp38d_n%d.mat',i);
    load(loadName);
    simRes1(i) = selResCell;
end
% simRes1
plotMode = 2; % multiple simulations
selResCell = simRes1;
referencePlotSimulations


%% let's then try to change the Keq to see if it could be
%
NumberCycles = 5;

% Start with change
% setup.clamp_GLCec = 0;
setup.clamp_GLCec = 1;
setup.changing_Keq_glt = 1;
tempValues = -[0 0.1 0.5 1 2 3];

% %% testing the change in the mass balance
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
% plotflag = 2; % variables by iteration
choosedataset

% loop to try different k_eq values
for i = 1:6
    %
    setup.Keq_glt_inc = tempValues(i);
    %
    setup.experiment = 1;
    [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    plotMode = 2; % multiple simulations
    % clear selResCell
    selResCell{1}.T_FF01 = T_FF01_1;
    selResCell{1}.Y_FF01 = Y_FF01_1;
    selResCell{1}.V_FF01 = V_FF01_1;
    referencePlotSimulations
end

%
plotMode = 0;
setup.clamp_GLCec = 0;






%% memoryDump

% %%
% % legendaMetabolites_addEnrichment;
% plotMode = 2; % multiple simulations
% referencePlotSimulations_enrichment
% plotMode = 0;


% %% checking sources of incoming labelled and non-labelled glucose
% temp_vGLT = selResCell{1}.V_FF01(:,1);
% temp_vNTH1 = selResCell{1}.V_FF01(:,20);
% temp_vATHv = selResCell{1}.V_FF01(:,46);
% temp_vGlycDeg = selResCell{1}.V_FF01(:,51);
% 
% temp_time = selResCell{1}.T_FF01;
% % glc_Ec (6 vs 65)
% temp_perc_labGLCec = selResCell{1}.Y_FF01(:,65) ./ selResCell{1}.Y_FF01(:,36);
% temp_perc_unlabGLCec = ones - temp_perc_labGLCec; 
% % tre_ic (25 vs 58)
% temp_perc_labTREic = selResCell{1}.Y_FF01(:,58) ./ selResCell{1}.Y_FF01(:,25);
% temp_perc_unlabTREic = ones - temp_perc_labTREic; 
% % tre_vac (38 vs 67)
% temp_perc_labTREvac = selResCell{1}.Y_FF01(:,67) ./ selResCell{1}.Y_FF01(:,38);
% temp_perc_unlabTREvac = ones - temp_perc_labTREvac; 
% % glyc_ic (40 vs 68)
% temp_perc_labGLYCic = selResCell{1}.Y_FF01(:,68) ./ selResCell{1}.Y_FF01(:,40);
% temp_perc_unlabGLYCic = ones - temp_perc_labGLYCic; 

% %%
% figure, plot(selResCell{1}.Y_FF01(:,6)), hold on
% plot(selResCell{1}.Y_FF01(:,65)),
% legend('6','65')
% %%
% figure,
% plot(temp_time,selResCell{1}.Y_FF01(:,6),'.-')
% hold on
% plot(temp_time,selResCell{1}.Y_FF01(:,65),'.-')
% legend('6','65')


% %%
% figure
% 
% subplot(2,6,1)
% plot(temp_time, temp_vGLT .* temp_perc_labGLCec,'.-')
% title('labeled via GLT')
% 
% subplot(2,6,2)
% plot(temp_time, temp_vNTH1 .* temp_perc_labTREic,'.-')
% title('labeled via NTH1')
% 
% subplot(2,6,3)
% plot(temp_time, temp_vATHv .* temp_perc_labTREvac,'.-')
% title('labeled via ATHv')
% 
% subplot(2,6,4)
% plot(temp_time, temp_vGlycDeg .* temp_perc_labGLYCic,'.-')
% title('labeled via GlycDeg')
% 
% subplot(2,6,6)
% temp1 = temp_vGLT .* temp_perc_labGLCec + temp_vNTH1 .* temp_perc_labTREic + temp_vATHv .* temp_perc_labTREvac + temp_vGlycDeg .* temp_perc_labGLYCic;
% plot(temp_time, temp1,'.-')
% title('labeled total')
% 
% 
% subplot(2,6,7)
% plot(temp_time,temp_vGLT .* temp_perc_unlabGLCec,'.-')
% title('unlabeled via GLT')
% 
% subplot(2,6,8)
% plot(temp_time, temp_vNTH1 .* temp_perc_unlabTREic,'.-')
% title('unlabeled via NTH1')
% 
% subplot(2,6,9)
% plot(temp_time, temp_vATHv .* temp_perc_unlabTREvac,'.-')
% title('unlabeled via ATHv')
% 
% subplot(2,6,10)
% plot(temp_time, temp_vGlycDeg .* temp_perc_unlabGLYCic,'.-')
% title('unlabeled via GlycDeg')
% 
% subplot(2,6,12)
% temp2 = temp_vGLT .* temp_perc_unlabGLCec + temp_vNTH1 .* temp_perc_unlabTREic + temp_vATHv .* temp_perc_unlabTREvac + temp_vGlycDeg .* temp_perc_unlabGLYCic;
% plot(temp_time, temp2,'.-')
% title('unlabeled total')
% 
% 
% suptitle('sources of labelled and unclabelled GLCic')


% %% correct
% figure,
% %
% subplot(1,2,1)
% plot(selResCell{1}.T_FF01, selResCell{1}.Y_FF01(:,26), '.-')
% title('t6p All')
% %
% subplot(1,2,2)
% plot(selResCell{1}.T_FF01, selResCell{1}.Y_FF01(:,59), '.-')
% title('t6p 13C labelled')
% suptitle('t6p, outside script')
% 
% %% here it's reversed (t6p-case)
% figure,
% %
% subplot(1,2,1)
% plot(T_FF01, Y_FF01(:,i), '.-')
% title('t6p 13C labelled')
% %
% subplot(1,2,2)
% plot(T_FF01, Y_FF01(:,correspMat(i)), '.-')
% title('t6p all')
% suptitle('t6p, in script')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %% simulation in the case of no enrichment
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% legendaMetabolites_addEnrichment;
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
% 
% %% tempPlotting -> based on initial simulations
% plotMode = 2;
% tempPlotting
% plotMode = 0;
% 
% % plotMode = 2;%11; % multiple simulations
% % referencePlotSimulations_enrichment
% % plotMode = 0;

