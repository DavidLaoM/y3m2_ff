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




%% Directly implementing the Csmin idea in out model
% general clamping and setup options
NumberCycles = 5;
setup.clamp_GLCec = 0;
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


% %% simulating enrichment
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


%% simulate with the latest paraemter estimates for PGI
% 
load('x102b.mat','xAll_102b');
xAll = [x; x; x; x; x];
xAll(2,57:60) = xAll_102b(2,57:60);
xAll(3,57:60) = xAll_102b(3,57:60);
xAll(4,57:60) = xAll_102b(4,57:60);
xAll(5,57:60) = xAll_102b(5,57:60);


%%  
legendaMetabolites_addEnrichment;
% plotflag = 2; % variables by iteration
choosedataset
for i = 1:5
    setup.experiment = 1;
    [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(xAll(i,:),canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    plotMode = 2; % multiple simulations
    % clear selResCell
    selResCell{1}.T_FF01 = T_FF01_1;
    selResCell{1}.Y_FF01 = Y_FF01_1;
    selResCell{1}.V_FF01 = V_FF01_1;
    % 
    referencePlotSimulations_enrichment
    plotMode = 0;
end



% setup.clamp_GLCec = 0;
% toc
% 
% 
% 
% % % % % %% simulating enrichment (modified)
% % % % % x2 = x; x2(164) = 0;
% % % % % setup.csmin = 0.094;%-0.05;
% % % % % 
% % % % % legendaMetabolites_addEnrichment;
% % % % % % plotflag = 2; % variables by iteration
% % % % % choosedataset
% % % % % setup.experiment = 1;
% % % % % % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % % plotMode = 2; % multiple simulations
% % % % % % clear selResCell
% % % % % selResCell{1}.T_FF01 = T_FF01_1;
% % % % % selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % % selResCell{1}.V_FF01 = V_FF01_1;
% % % % % % 
% % % % % referencePlotSimulations_enrichment
% 
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %% Parameter estimation setup
% 
% % blank and constant setup
% legendaMetabolites_addEnrichment;
% choosedataset
% setup.ExpData = ExpData;
% % 
% blankWeight = zeros(1,89); % 85+3 for glycerol, 89 for glucose enrichment
% lambdalist = 0;
% setup.parEst.lambda = lambdalist(1); lam = setup.parEst.lambda;
% % % % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx});
% % % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
% % %     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% % options = optimoptions('lsqnonlin','Display','iter',...
% %     'OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
% %     'FiniteDifferenceStepSize',0.2);
% options = optimoptions('lsqnonlin','Display','iter',...
%     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4,...
%     'DiffMinChange',0.1);
% %     'OutputFcn',{@saveIterationsMain},...
% %     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% % setup.parEst.costfun = 10015; 
% 
% % % optimOptions2
% % enzNum = 6;
% % nchoosek_cell = cell(enzNum,1);
% % for i = 1:enzNum
% %     nchoosek_cell{i} = nchoosek(1:enzNum,i);
% % end
% % % concatenate
% % if exist('z','var')
% %     clear z
% % end
% % for i = 1:(enzNum-1)
% %     if exist('z','var')
% %     else
% %         z = nchoosek_cell{i}';
% %     end
% %     B = nchoosek_cell{i+1}';
% %     sA = size(z);
% %     sB = size(B);
% %     C = max(sA(1),sB(1));
% %     z = [[z;zeros(abs([C 0]-sA))],[B;zeros(abs([C,0]-sB))]];
% % end
% % z = z';
% % % z = z([1:3,5:6,8,11,15],:);
% 
% % % nParCombs = 1;
% % [nParCombs,~] = size(z);
% % nWeightCombs = 8; % 3 specific + 5 growing part
% % nExtras = 0;
% % nMS = 0;
% % ntests = nParCombs * nWeightCombs + nExtras;
% ntests = 15;
% nMS = 20;
% 
% % % names cell (to save)
% % names_cell = cell(ntests,1);
% % tempName = 'FF_pE2_pCs%d_wA%d'; %.mat to be added upon call
% % for i = 1:ntests
% %     pC_num = fix((i-1)/nWeightCombs) + 1;
% %     warr_num = rem(i,nWeightCombs);
% %     if warr_num == 0, warr_num = nWeightCombs; end
% %     names_cell{i} = sprintf(tempName,pC_num,warr_num);
% % end
% % optimOptions2.names_cell = names_cell;
% names_cell = cell(ntests,1);
% % for i = 1:ntests
% %     names_cell{i} = sprintf('FF_pE9_wA%d',i);
% % end
% % 
% names_cell{1} = 'FF_pE13b_pC1_wA1_nS%d.mat';
% names_cell{2} = 'FF_pE13b_pC1_wA2_nS%d.mat';
% names_cell{3} = 'FF_pE13b_pC1_wA3_nS%d.mat';
% names_cell{4} = 'FF_pE13b_pC1_wA4_nS%d.mat';
% names_cell{5} = 'FF_pE13b_pC1_wA5_nS%d.mat';
% % 
% names_cell{6} = 'FF_pE13b_pC2_wA1_nS%d.mat';
% names_cell{7} = 'FF_pE13b_pC2_wA2_nS%d.mat';
% names_cell{8} = 'FF_pE13b_pC2_wA3_nS%d.mat';
% names_cell{9} = 'FF_pE13b_pC2_wA4_nS%d.mat';
% names_cell{10} = 'FF_pE13b_pC2_wA5_nS%d.mat';
% % 
% names_cell{11} = 'FF_pE13b_pC3_wA1_nS%d.mat';
% names_cell{12} = 'FF_pE13b_pC3_wA2_nS%d.mat';
% names_cell{13} = 'FF_pE13b_pC3_wA3_nS%d.mat';
% names_cell{14} = 'FF_pE13b_pC3_wA4_nS%d.mat';
% names_cell{15} = 'FF_pE13b_pC3_wA5_nS%d.mat';
% %
% %
% optimOptions2.names_cell = names_cell;
% 
% % par combinations cell: to select the parameters to optimize
% %
% parsGLT = [35 36 38]; % glt
% parsHXK = 28:34; % hxk
% %
% parsPGM1 = 83:86; % pgm1
% parsUGP = 144:148; % ugp
% parsTPS1 = 124:128; % tps1
% parsTPS2 = 119:121; % tps2
% parsNTH1 = 122:123; % nth1
%     parsTreCycle = [parsPGM1, parsUGP, parsTPS1, parsTPS2, parsNTH1];
% % parsGLY = 159;
% parsGLY2 = 160:161;
%     %
%     parsATH1 = [149:150,162:165];       % ath1_v.e  [149:150]
%         parsATH1vac = [149 150 162];
%         parsATH1ec = [163 164 165];
% %     parsAGT1 = [151:154,158,166];   % agt1      [151:154,158]
%     parsAGT1 = [151:154,158];   % agt1      [151:154,158]
%     parsVACT = 155:156;       % vacT      [155:156]
%     parsTreExport = [parsATH1, parsAGT1, parsVACT];
% 
% pC1 = [parsGLT, parsATH1ec];
% pC2 = [parsGLT, parsHXK, parsTreExport];
% pC3 = [parsGLT, parsHXK, parsTreCycle, parsTreExport];
% %
% parComb_cell = cell(ntests,1);
% % 
% parComb_cell{1} = pC1;
% parComb_cell{2} = pC1;
% parComb_cell{3} = pC1;
% parComb_cell{4} = pC1;
% parComb_cell{5} = pC1;
% % 
% parComb_cell{6} = pC2;
% parComb_cell{7} = pC2;
% parComb_cell{8} = pC2;
% parComb_cell{9} = pC2;
% parComb_cell{10} = pC2;
% % 
% parComb_cell{11} = pC3;
% parComb_cell{12} = pC3;
% parComb_cell{13} = pC3;
% parComb_cell{14} = pC3;
% parComb_cell{15} = pC3;
% % 
% 
% % weights:
% w_glc_ic = 6;
% w_glc_ec = 36;
% w_g1p = 21;
% w_udpg = 24;
% w_t6p_tre = [26 25];
% w_treRates_glk_pgi = [55 59 57 58 56 40 41];
% w_glt = 39;
% w_gly = [86 87 88];
% w_tre_ic = 25;
% w_tre_ec = 37;
% w_ath1 = [83 84];
%     w_ath1_ec = 83;
%     w_ath1_vac = 84; 
% w_agt1 = 85;
% w_glc_enrich = 89;
% %
% w_idxs1 = w_glc_enrich;
% w_idxs2 = [w_glc_enrich, w_glc_ec];
% w_idxs3 = [w_glc_enrich, w_glt];
% w_idxs4 = [w_glc_enrich, w_glc_ec, w_glt];
% w_idxs5 = [w_glc_enrich, w_glc_ec, w_glt, ...
%     w_treRates_glk_pgi,...
%     w_g1p, w_udpg, w_t6p_tre, w_tre_ic, w_tre_ec];
% %
% warray_cell = cell(ntests,1);
% for i = 1:ntests % base,
%     warray_cell{i} = blankWeight;
% end
% % 
% warray_cell{1}(w_idxs1) = ones;
% warray_cell{2}(w_idxs2) = ones;
% warray_cell{3}(w_idxs3) = ones;
% warray_cell{4}(w_idxs4) = ones;
% warray_cell{5}(w_idxs5) = ones;
% % 
% warray_cell{6}(w_idxs1) = ones;
% warray_cell{7}(w_idxs2) = ones;
% warray_cell{8}(w_idxs3) = ones;
% warray_cell{9}(w_idxs4) = ones;
% warray_cell{10}(w_idxs5) = ones;
% % 
% warray_cell{11}(w_idxs1) = ones;
% warray_cell{12}(w_idxs2) = ones;
% warray_cell{13}(w_idxs3) = ones;
% warray_cell{14}(w_idxs4) = ones;
% warray_cell{15}(w_idxs5) = ones;
% % 
% 
% 
% %% Selected runs
% % 1*165
% selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS
% selectedRuns_idx = selCase
% % if selCase <= 24
% %     selectedRuns_idx = [1,17,33] + 1 * (selCase - 1)
% % else
% %     selectedRuns_idx = [1:15] + 15 * (selCase - 25) + 384;
% % end
% % disp(selectedRuns_idx)
% % % % % selectedRuns_idx = 1:ntests;
% % 
% % 
% %% PARAMETER ESTIMATION: NO MS
% for o = selectedRuns_idx
%     
%     for u = 1:nMS
%         % disp
%         fprintf('o = %f, u = %f.\n', o, u) %disp(o)
% 
%         % find in the precreated arrays
%         selPars = parComb_cell{o};
%         warray = warray_cell{o}; 
%         tempName = names_cell{o};
%         saveName = sprintf(tempName,u);
% 
%         % core options:
%         x_iter = x;
%         if u == 1 % initial
%         elseif u == 2 % ath1ec hydrolysis
%             x_iter(164) = 0;
%         end
%         x_temp = x_iter(selPars);
%         if u > 2 % rand
%             rng(u); x_temp = x_temp - 1 + 2 * rand(size(x_temp));%x_iter
%         end
%         plength = length(selPars);
%         lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
%         ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;
%         if(o >= 6)
%             lb(end-12:end) = -10;
%             ub(end-12:end) = 10;
%         else
%             lb(end-2:end) = -10;
%             ub(end-2:end) = 10;
%         end
% 
%     %     % %% run check
%     %     [error]=costfunSystemY3M1_FF_pRun(x_temp,canelas_SS,setup,xarray,data,dataset,NumberCycles,IC,selPars,warray);
%     %     errorAnalysis_Y3M1;
%     %     % %%
% 
% %         %
% %         [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x_iter,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% %         plotMode = 2; % multiple simulations
% %         % clear selResCell
% %         selResCell{1}.T_FF01 = T_FF01_1;
% %         selResCell{1}.Y_FF01 = Y_FF01_1;
% %         selResCell{1}.V_FF01 = V_FF01_1;
% %         % 
% %         referencePlotSimulations_enrichment
% %         plotMode = 0;
% %         %
%         
%         % parameter estimation
%         sprintf('run num.%d tested',o)
%         tic
% %         [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun,x_temp,lb,ub,options,canelas_SS,setup,x_iter,data,dataset,NumberCycles,IC,selPars,warray);
%         [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun_enrichment,x_temp,lb,ub,options,canelas_SS,setup,x_iter,data,dataset,NumberCycles,IC,selPars,warray);
%         t = toc; 
%         disp(xres);
% 
%         % saving
%     % % % %     parsave_Y3M2_cluster(saveName, xres, resnorm, residual, exitflag, warray, t, selPars);
%         save(saveName, 'xres', 'resnorm', 'residual', 'exitflag', 'warray', 't', 'selPars');
%     end
%     
% end
% % disp('parallel run below')
% % % save
% % saveName = 'FF_pE4_x1.mat';
% % save(saveName,'allRes_xres','allRes_resnorm','allRes_residual','allRes_exitflag','allRes_t','allRes_warray','allRes_selPars');
% % % % % delete(pool)
% % quit
% 
% 
% % % % % %% DATA RECALL (not the multistart cases yet)
% % % % % % xFullScale_99b_parCombs5_warray8.mat is lacking at start
% % % % % xAll1 = x;
% % % % % namesHits1 = cell(1,1); namesHits1{1} = 'initial';
% % % % % for i = 1:ntests
% % % % %     loadName = names_cell{i};
% % % % %     if exist(loadName,'file') == 2 % if exist
% % % % % % % % %             loadName = [folder, '\workContinuation\tempResults\', loadName];
% % % % %         load(loadName);
% % % % %         selPars = parComb_cell{i};
% % % % %         x3 = x;
% % % % %         x3(selPars) = xres;
% % % % %         xAll1 = [xAll1; x3];
% % % % %         namesHits1 = [namesHits1; loadName];
% % % % % % % % %             namesHits1 = [namesHits1; loadName(end-21:end)];
% % % % %         
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
% % % % % %% simulate results
% % % % %     % 1-6
% % % % %     % 7-14
% % % % % nParts = length(namesHits1);
% % % % % simRes1 = cell(1,nParts);
% % % % % parpool(4)
% % % % % parfor i = 1:nParts
% % % % %     disp(i);
% % % % %     xSel = xAll1(i,:);
% % % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % %     simRes1{i}.T_FF01 = T_FF01;
% % % % %     simRes1{i}.Y_FF01 = Y_FF01;
% % % % %     simRes1{i}.V_FF01 = V_FF01;
% % % % % end
% % % % % save('pE11f_simRes.mat','simRes1')
% % % % % % % quit
% % % % % % % save('pE8_simRes_updated.mat','simRes1')
% % % % % %%
% % % % % nParts = length(namesHits1);
% % % % % simRes1 = cell(1,nParts);
% % % % % load('pE11f_simRes.mat','simRes1')
% % % % % 
% % % % % 
% % % % % %% plot data
% % % % % [legenda] = legendaFull; %legenda for the names needed
% % % % % metNames = legenda.metabolites;
% % % % % reactNames = legenda.fluxes;
% % % % % 
% % % % % choosedataset
% % % % % 
% % % % % % plotMode = 1; % single simulation
% % % % % plotMode = 2; % multiple simulations
% % % % % % plotMode = 10; % all
% % % % % 
% % % % % % tempSave for the entire thing
% % % % % tempSave = simRes1;
% % % % % 
% % % % % % pC1
% % % % % selResCell = tempSave(1:4);
% % % % % referencePlotSimulations
% % % % % 
% % % % % % pC2
% % % % % selResCell = tempSave([1,5:7]);
% % % % % referencePlotSimulations
% % % % % 
% % % % % % pC3
% % % % % selResCell = tempSave([1,8:10]);
% % % % % referencePlotSimulations
% % % % % 
% % % % % % pC4
% % % % % selResCell = tempSave([1,11:13]);
% % % % % referencePlotSimulations
% % % % % 
% % % % % % pC5
% % % % % selResCell = tempSave([1,14:16]);
% % % % % referencePlotSimulations
% % % % % 
% % % % % % pC6
% % % % % selResCell = tempSave([1,17]);
% % % % % referencePlotSimulations
% % % % % 
% % % % % % % % all sims
% % % % % % % % simRes1 = tempSave;
% % % % % % % selResCell = simRes1;
% % % % % % % referencePlotSimulations
% % % % % 
% % % % % % %%
% % % % % for i = 1:15
% % % % %     fid = fopen('parEst13b_GLTcsmin_parEst.m','rt');
% % % % %     X = fread(fid);
% % % % %     fclose(fid);
% % % % %     X = char(X.');
% % % % % %     str2rep = sprintf('vals2run = %d; % <--',i);
% % % % %     str2rep = sprintf('selCase = %d; % % <== changed already',i);
% % % % %     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
% % % % %     tempName = sprintf('pE13b_%d.m',i);
% % % % %     fid2 = fopen(tempName,'wt') ;
% % % % %     fwrite(fid2,Y) ;
% % % % %     fclose(fid2);
% % % % % end



