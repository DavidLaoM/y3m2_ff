% % sim_pEST_FF.m
% This is the main file to make changes on the model, estimate parameters
% and simulate to see the output / study the model.
% This file is organised in the following sections:
% (0) Initial setup
    % The whole model is quite complex and needs multiplie inputs to run. In
    % this project, we'll leave most constant and just focus in the regions of
    % interest. This section is to be left as it is, for now.
% (1) simulations (GP suarez)
    % Simulations of the 3 datasets to be used in this study. To determine
    % which dataset to work with, the combinations have to be used:
        % Dataset 1: setup.GPdataset.GP400WT = 1; setup.GPdataset.GP1800WT = 0; setup.GPdataset.GP400M = 0;
        % Dataset 2: setup.GPdataset.GP400WT = 0; setup.GPdataset.GP1800WT = 1; setup.GPdataset.GP400M = 0;
        % Dataset 3: setup.GPdataset.GP400WT = 0; setup.GPdataset.GP1800WT = 0; setup.GPdataset.GP400M = 1;
% (2) Parameter estimation and resulting model simulations
    % Here is the area of study. You can use the legend to find the names
    % of what you want to study.
    % First the dataset to study can be selected: (D1) 1-0-0, (D3) 0-1-0 or (D4) 0-0-1
        % % select dataaset to optimize for:
        % setup.GPdataset.GP400WT = 1;
        % setup.GPdataset.GP1800WT = 0;
        % setup.GPdataset.GP400M = 0;
    % Select the paraemters of study by their index. By default constrained
    % to glt and hxk. The longer the array of parameters selected (the more
    % parameters that are taken into account at the same time) the higher
    % the computational cost (slower to give results).
        % % select parameters to study:
        % setup.caseStudy.parameters = [(28:34),(35:36)]; % HXK: 28:34, GLT: 35:36.
    % Select what you consider in the cost function. You always have to
    % tell the algorithm what to select to optimize for. In this case, the
    % array 'w' has a component for each metabolite, to be selected by
    % index. 
    % For example, if all is zeros, no optimization happens, but if w(5) is
    % 1, it will optimize for G6P (metabolite #5 in the legend). If w(4) =
    % 2 and w(5) = 1, it will optimize for both F6P and G6P, but give
    % double importance to F6P than G6P.
        % % select weights to optimize:
        % w = zeros(32,1); w(5) = 0;
        % setup.w = w;
    % To use later. Regularization factor. If the absolute value increases,
    % the algorithm will push so that parameters are close to the
    % literature value (=0).
        % % select regularization factor:
        % setup.parEst.lambda = 0;

%% (0) Initial setup.

% startup
clear, close all
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

load('x_Comb32_flx');

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

%% updated parameter sets
load('x_Comb32_flx'); x = xarray;
load('pset_20210325.mat')
load('parSet_99b.mat')
%
parsGlyco = [28:38];
parsTre = [83:86,119:128];
parsAXP = [109:111,129:131];
parsAdd = 144:158;
%
x1 = x;
x2 = x99b;
    x2(parsGlyco) = x(parsGlyco);
    x2(parsTre) = x(parsTre);
%     x2(parsAXP) = x(parsAXP);
    x2(parsAdd) = x(parsAdd);
x3 = x2;
    x3([109 129]) = x([109 129]);



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
% plotflag = 2; % variables by iteration
choosedataset
setup.experiment = 1;
% [T_FF01,Y_FF01,V_FF01] = simulate_FF(x3,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);


%% Check that other transport rates are well set to zero
parsFRT = 137:139;
parsHXK_FRT = 140:143;


%% PART 1. SCATTER OVER ALL PARAMETERS
nS = 100;
optionSample = 1; % % % % Change here for parallel
% randArray is a big array from which error subsection is later selected
rng(optionSample); randArray = lhsdesign(nS,158) * 2 - 1;


%% SAMPLING: samples creation
% pars range [-0.5 +0.5]
parsGLT = [35 36 38];
parsGLK = 28:34;
% pars range [-1.5 +1.5]
parsPGM1 = 83:86;
parsTPS1 = 124:128;
parsTPS2 = 119:121;
parsNTH1 = 122:123;
parsUGP = 144:148;
% pars range [-3 +3]
parsATH1 = 149:150;
parsAGT1 = [151:154,158];
parsvacT = 155:157;

% creating temp_xArray
xScatter = zeros(nS, length(x3));
for i = 1:nS
    xScatter(i,:) = x3;
end
xScatter(:,[parsGLT parsGLK]) = xScatter(:,[parsGLT parsGLK]) + 0.5*randArray(:,[parsGLT parsGLK]);
xScatter(:,[parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP]) = xScatter(:,[parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP]) + 1.5*randArray(:,[parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP]);
xScatter(:,[parsATH1 parsAGT1 parsvacT]) = xScatter(:,[parsATH1 parsAGT1 parsvacT]) + 3.0*randArray(:,[parsATH1 parsAGT1 parsvacT]);

% selection and editing
% some enzymes
if optionSample == 1, selPars = [parsGLT parsGLK]; %end
elseif optionSample == 2, selPars = parsGLT; %end
elseif optionSample == 3, selPars = parsGLK; %end
%
elseif optionSample == 4, selPars = [parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP]; %end
elseif optionSample == 5, selPars = parsPGM1; %end
elseif optionSample == 6, selPars = parsTPS1; %end
elseif optionSample == 7, selPars = parsTPS2; %end
elseif optionSample == 8, selPars = parsNTH1; %end
elseif optionSample == 9, selPars = parsUGP; %end
%
elseif optionSample == 10, selPars = [parsATH1 parsAGT1 parsvacT]; %end
elseif optionSample == 11, selPars = parsATH1; %end
elseif optionSample == 12, selPars = parsAGT1; %end
elseif optionSample == 13, selPars = parsvacT; %end
% all enZymes
else, selPars = [parsGLT parsGLK parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP parsATH1 parsAGT1 parsvacT]; end


%% SAMPLING: simulations and storage
% saveName_1 = sprintf('saveRes_section1_option%d.mat',optionSample);
% saveName_2 = sprintf('saveRes_section2_option%d.mat',optionSample);
% saveName_3 = sprintf('saveRes_section3_option%d.mat',optionSample);
% saveName_4 = sprintf('saveRes_section4_option%d.mat',optionSample);
% saveName_5 = sprintf('saveRes_section5_option%d.mat',optionSample);
% 
% nParts = 20;
% % part 1
% tempSave1 = cell(1,nParts);
% for i = 1:nParts
%     i2 = i + 0;
%     xSel = x3;
%     xSel(selPars) = xScatter(i2,selPars);
% %     T_FF01 = 'dummyTsave';
% %     Y_FF01 = 'dummyYsave';
% %     V_FF01 = 'dummyVsave';
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     tempSave1{i}.T_FF01 = T_FF01;
%     tempSave1{i}.Y_FF01 = Y_FF01;
%     tempSave1{i}.V_FF01 = V_FF01;
% end
% save(saveName_1,'tempSave1'); clear tempSave1
% % part 2
% tempSave2 = cell(1,nParts);
% for i = 1:nParts
%     i2 = i + 20;
%     xSel = x3;
%     xSel(selPars) = xScatter(i2,selPars);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     tempSave2{i}.T_FF01 = T_FF01;
%     tempSave2{i}.Y_FF01 = Y_FF01;
%     tempSave2{i}.V_FF01 = V_FF01;
% end
% save(saveName_2,'tempSave2'); clear tempSave2
% % part 3
% tempSave3 = cell(1,nParts);
% for i = 1:nParts
%     i2 = i + 40;
%     xSel = x3;
%     xSel(selPars) = xScatter(i2,selPars);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     tempSave3{i}.T_FF01 = T_FF01;
%     tempSave3{i}.Y_FF01 = Y_FF01;
%     tempSave3{i}.V_FF01 = V_FF01;
% end
% save(saveName_3,'tempSave3'); clear tempSave3
% % part 4
% tempSave4 = cell(1,nParts);
% for i = 1:nParts
%     i2 = i + 60;
%     xSel = x3;
%     xSel(selPars) = xScatter(i2,selPars);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     tempSave4{i}.T_FF01 = T_FF01;
%     tempSave4{i}.Y_FF01 = Y_FF01;
%     tempSave4{i}.V_FF01 = V_FF01;
% end
% save(saveName_4,'tempSave4'); clear tempSave4
% % part 5
% tempSave5 = cell(1,nParts);
% for i = 1:nParts
%     i2 = i + 80;
%     xSel = x3;
%     xSel(selPars) = xScatter(i2,selPars);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     tempSave5{i}.T_FF01 = T_FF01;
%     tempSave5{i}.Y_FF01 = Y_FF01;
%     tempSave5{i}.V_FF01 = V_FF01;
% end
% save(saveName_5,'tempSave5'); clear tempSave5


%% recall data
selPars_cell = cell(1,23);
for i = 1:23
    % some enzymes
    if i == 1, selPars_cell{i} = [parsGLT parsGLK]; %end
    elseif i == 2, selPars_cell{i} = parsGLT; %end
    elseif i == 3, selPars_cell{i} = parsGLK; %end
    %
    elseif i == 4, selPars_cell{i} = [parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP]; %end
    elseif i == 5, selPars_cell{i} = parsPGM1; %end
    elseif i == 6, selPars_cell{i} = parsTPS1; %end
    elseif i == 7, selPars_cell{i} = parsTPS2; %end
    elseif i == 8, selPars_cell{i} = parsNTH1; %end
    elseif i == 9, selPars_cell{i} = parsUGP; %end
    %
    elseif i == 10, selPars_cell{i} = [parsATH1 parsAGT1 parsvacT]; %end
    elseif i == 11, selPars_cell{i} = parsATH1; %end
    elseif i == 12, selPars_cell{i} = parsAGT1; %end
    elseif i == 13, selPars_cell{i} = parsvacT; %end
    % all enZymes
    else, selPars_cell{i} = [parsGLT parsGLK parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP parsATH1 parsAGT1 parsvacT]; end
end

fullSaveRes = cell(23,100);
for i = 1:23
    % names to load
    loadName_1 = sprintf('saveRes_section1_option%d.mat',i);
    loadName_2 = sprintf('saveRes_section2_option%d.mat',i);
    loadName_3 = sprintf('saveRes_section3_option%d.mat',i);
    loadName_4 = sprintf('saveRes_section4_option%d.mat',i);
    loadName_5 = sprintf('saveRes_section5_option%d.mat',i);
    % load and place
    load(loadName_1);
    tempLoadRes = tempSave1;
    load(loadName_2);
    tempLoadRes = [tempLoadRes, tempSave2];
    load(loadName_3);
    tempLoadRes = [tempLoadRes, tempSave3];
    load(loadName_4);
    tempLoadRes = [tempLoadRes, tempSave4];
    load(loadName_5);
    tempLoadRes = [tempLoadRes, tempSave5];
    % saveInFullLocation
    fullSaveRes(i,:) = tempLoadRes;    
end

% %%
% nlen = 100;
% len_overall = zeros(1,nlen);
% max_overall = zeros(1,nlen);
% for i = 1:nlen
%     len_overall(i) = length(fullSaveRes{1,i}.T_FF01);
%     max_overall(i) = max(fullSaveRes{1,i}.T_FF01);
% end
% %%
% figure,
% subplot(1,2,1), plot(len_overall,'.-')
% subplot(1,2,2), plot(max_overall,'.-')

%% visualization scatter simulations
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;

choosedataset

% plotMode = 1; % single simulation
plotMode = 2; % multiple simulations
% plotMode = 10; % all
for o = 3 % selPars = [parsGLT parsGLK]; %end
% for o = 4:9 % selPars = [parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP]; %end
% for o = 10:13 % selPars = [parsATH1 parsAGT1 parsvacT]; %end
% for o = 14:23 % all
    selResCell = fullSaveRes(o,:);
    referencePlotSimulations
end


%% Correlation: selResCell
% recall
pep_num8_min = zeros(nS,1);
v_glk_dp3 = zeros(nS,1);
v_glk_dp4 = zeros(nS,1);
for i = 1:nS
    pep_num8_min(i) = interp1(selResCell{i}.T_FF01,selResCell{i}.Y_FF01(:,12),60,'pchip');
    v_glk_dp3(i) = interp1(selResCell{i}.T_FF01,selResCell{i}.V_FF01(:,2),80,'pchip');
    v_glk_dp4(i) = interp1(selResCell{i}.T_FF01,selResCell{i}.V_FF01(:,3),220,'pchip');
end
% plotting
figure
subplot(2,1,1)
scatter(v_glk_dp3,pep_num8_min,5,'black','filled')
xlabel('glk_{3}'), ylabel('pep_8')
subplot(2,1,2)
scatter(v_glk_dp4,pep_num8_min,5,'black','filled')
xlabel('glk_{4}'), ylabel('pep_8')


%% Correlation: reorder
[~,idxs_reorder] = sort(v_glk_dp4);

selResCell = fullSaveRes(3,idxs_reorder);
referencePlotSimulations

%%
[T_FF01,Y_FF01,V_FF01] = simulate_FF(x3,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
clear selResCell
selResCell{1}.T_FF01 = T_FF01;
selResCell{1}.Y_FF01 = Y_FF01;
selResCell{1}.V_FF01 = V_FF01;
%     selResCell = fullSaveRes(o,:);
referencePlotSimulations
% end

%% memorySave
% 
% 
% %% PART 2 SIMULATIONS WITH TREHALOSE METABOLISM TURNED OFF
% 
% xSel_start = x3;
% xSel_noTreComp = x3;
%     xSel_noTreComp([510 152 155]) = -10 * ones;
% 
% tempSaveSilenceTreCompart = cell(1,2);
% % 1
% [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel_start,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% tempSaveSilenceTreCompart{1}.T_FF01 = T_FF01;
% tempSaveSilenceTreCompart{1}.Y_FF01 = Y_FF01;
% tempSaveSilenceTreCompart{1}.V_FF01 = V_FF01;
% % 2
% [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel_noTreComp,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% tempSaveSilenceTreCompart{2}.T_FF01 = T_FF01;
% tempSaveSilenceTreCompart{2}.Y_FF01 = Y_FF01;
% tempSaveSilenceTreCompart{2}.V_FF01 = V_FF01;
% 
% %%
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% 
% choosedataset
% 
% % plotMode = 1; % single simulation
% plotMode = 2; % multiple simulations
% % plotMode = 10; % all
% % % % % for o = 1:2
%     selResCell = tempSaveSilenceTreCompart(1,:);
%     referencePlotSimulations
% % % % % end
% % end
% 
% 
% %% PART 3. ANALYZING THE SCATTER PLOT FOR THE 'GLYCOGEN ZONE'
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% 
% choosedataset
% 
% % plotMode = 1; % single simulation
% plotMode = 2; % multiple simulations
% % plotMode = 10; % all
% % for o = 1:3 % selPars = [parsGLT parsGLK]; %end
% for o = 4:9 % selPars = [parsPGM1 parsTPS1 parsTPS2 parsNTH1 parsUGP]; %end
% % for o = 10:13 % selPars = [parsATH1 parsAGT1 parsvacT]; %end
% % for o = 14:23 % all
%     selResCell = fullSaveRes(o,:);
%     referencePlotSimulations
% end
% 
% 
% %% correlation plots
% corrrelation_plots_setup.numPlots = 2;
% corrrelation_plots_setup.selData = cell(1,corrrelation_plots_setup.numPlots);
%     corrrelation_plots_setup.selData{1} = [21 24 26 25];
%     corrrelation_plots_setup.selData{2} = [17 18 21 19 20];
% corrrelation_plots_setup.refFlux = 3;
% corrrelation_plots_setup.criteria = 'maxValue_normVals';
% 
% corrrelation_plots_setup.numGroups = 5;
% corrrelation_plots_setup.groups = [5 6 7 8 9];
% % corrrelation_plots_setup.groups = [4 5 6 7 8 9]; 
% 
% % corrrelation_plots_setup.labelsGroups = [5 6 7 8 9]; 
% corrrelation_plots_setup.labelsGroups = {...
%     'G5, pgm1',...
%     'G6, tps1',...
%     'G7, tps2',...
%     'G8, nth1',...
%     'G9, ugp'};
% corrrelation_plots_setup.labels_selData1 = {...
%     'vPGI',...
%     'g1p',...
%     'udpg',...
%     't6p',...
%     'tre_{cyt}'};
% corrrelation_plots_setup.labels_selData2 = {...
%     'vPGI',...
%     'pgm1',...
%     'ugp',...
%     'tps1',...
%     'tps2',...
%     'nth1'};
% 
% 
% %%
% assay1001_correlationInSimulations;
% % 

%%
% % % %%
% for i = 1:23
%     fid = fopen('assay1_scatter_initial_state.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('optionSample = %d % % <== changed already',i);
%     Y = strrep(X,'optionSample = 1; % % % % Change here for parallel',str2rep);
%     tempName = sprintf('ff_scatter1_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end


