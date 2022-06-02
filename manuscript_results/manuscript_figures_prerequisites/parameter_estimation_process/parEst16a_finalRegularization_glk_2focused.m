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

%% (2021 - 08 - 09) right balance around ATH1ec
setup.TRE_recirculation_rightBalances = 1;


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
choosedataset


%% Starting simulation
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % 
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;
% % 
% referencePlotSimulations_enrichment
% plotMode = 0;


%% Recall parameter sets
% safecopy the simulations
safecopy_startSim = selResCell{1};
% Current parameter set 'x_Y3M2' and previous to regularize for 'x_Y3M1'.
load('parSet_99b.mat'); x_Y3M1 = x99b;
x_Y3M2 = x;
% region of interest
parsGLT = 35:38;
parsGLK = 28:34;
parsTre = [83:86,119:128];
parsAXP = [109:111,129:131];


%% Parameter estimation setup
% blank and constant setup
legendaMetabolites_addEnrichment;
choosedataset
setup.ExpData = ExpData;
setup.simData = safecopy_startSim;
setup.pset_Y3M1 = x_Y3M1;
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

% % optimOptions2
% enzNum = 6;
% nchoosek_cell = cell(enzNum,1);
% for i = 1:enzNum
%     nchoosek_cell{i} = nchoosek(1:enzNum,i);
% end
% % concatenate
% if exist('z','var')
%     clear z
% end
% for i = 1:(enzNum-1)
%     if exist('z','var')
%     else
%         z = nchoosek_cell{i}';
%     end
%     B = nchoosek_cell{i+1}';
%     sA = size(z);
%     sB = size(B);
%     C = max(sA(1),sB(1));
%     z = [[z;zeros(abs([C 0]-sA))],[B;zeros(abs([C,0]-sB))]];
% end
% z = z';
% % z = z([1:3,5:6,8,11,15],:);

% % nParCombs = 1;
% [nParCombs,~] = size(z);
% nWeightCombs = 8; % 3 specific + 5 growing part
% nExtras = 0;
% nMS = 0;
% ntests = nParCombs * nWeightCombs + nExtras;
ntests = 16;

% % names cell (to save)
% names_cell = cell(ntests,1);
% tempName = 'FF_pE2_pCs%d_wA%d'; %.mat to be added upon call
% for i = 1:ntests
%     pC_num = fix((i-1)/nWeightCombs) + 1;
%     warr_num = rem(i,nWeightCombs);
%     if warr_num == 0, warr_num = nWeightCombs; end
%     names_cell{i} = sprintf(tempName,pC_num,warr_num);
% end
% optimOptions2.names_cell = names_cell;
names_cell = cell(ntests,1);
% for i = 1:ntests
%     names_cell{i} = sprintf('FF_pE9_wA%d',i);
% end
% 
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
% 
optimOptions2.names_cell = names_cell;

% par combinations cell: to select the parameters to optimize
%
parsGLT = [35 36 38]; % glt
parsHXK = 28:34; % hxk
%
parsPGM1 = 83:86; % pgm1
parsUGP = 144:148; % ugp
parsTPS1 = 124:128; % tps1
parsTPS2 = 119:121; % tps2
parsNTH1 = 122:123; % nth1
    parsTreCycle = [parsPGM1, parsUGP, parsTPS1, parsTPS2, parsNTH1];
% parsGLY = 159;
parsGLY2 = 160:161;
    %
    parsATH1 = [149:150,162:165];       % ath1_v.e  [149:150]
        parsATH1vac = [149 150 162];
        parsATH1ec = [163 164 165];
%     parsAGT1 = [151:154,158,166];   % agt1      [151:154,158]
    parsAGT1 = [151:154,158];   % agt1      [151:154,158]
    parsVACT = 155:156;       % vacT      [155:156]
    parsTreExport = [parsATH1, parsAGT1, parsVACT];

pC16a = parsHXK;
pC16b = [parsGLT, parsHXK, parsTreCycle];
%
parComb_cell = cell(ntests,1);
% 
parComb_cell{1} = pC16a;
parComb_cell{2} = pC16a;
parComb_cell{3} = pC16a;
parComb_cell{4} = pC16a;
parComb_cell{5} = pC16a;
parComb_cell{6} = pC16a;
parComb_cell{7} = pC16a;
parComb_cell{8} = pC16a;
parComb_cell{9} = pC16a;
parComb_cell{10} = pC16a;
parComb_cell{11} = pC16a;
parComb_cell{12} = pC16a;
parComb_cell{13} = pC16a;
parComb_cell{14} = pC16a;
parComb_cell{15} = pC16a;
parComb_cell{16} = pC16a;
% 
% % % % parComb_cell{1} = pC16b;
% % % % parComb_cell{2} = pC16b;
% % % % parComb_cell{3} = pC16b;
% % % % parComb_cell{4} = pC16b;
% % % % parComb_cell{5} = pC16b;
% % % % parComb_cell{6} = pC16b;
% % % % parComb_cell{7} = pC16b;
% % % % parComb_cell{8} = pC16b;
% % % % parComb_cell{9} = pC16b;
% % % % parComb_cell{10} = pC16b;
% % % % parComb_cell{11} = pC16b;
% % % % parComb_cell{12} = pC16b;
% % % % parComb_cell{13} = pC16b;
% % % % parComb_cell{14} = pC16b;
% % % % parComb_cell{15} = pC16b;
% % % % parComb_cell{16} = pC16b;
% 

% weights: 
% + w(60) -> simulated hxk.
% + w(61) -> reference parameter set.
w_simHXK = 60;
w_refY3M1 = 61;
% % w_glc_ic = 6;
% % w_glc_ec = 36;
% % w_g1p = 21;
% % w_udpg = 24;
% % w_t6p_tre = [26 25];
% % w_treRates_glk_pgi = [55 59 57 58 56 40 41];
% % w_glt = 39;
% % w_gly = [86 87 88];
% % w_tre_ic = 25;
% % w_tre_ec = 37;
% % w_ath1 = [83 84];
% %     w_ath1_ec = 83;
% %     w_ath1_vac = 84; 
% % w_agt1 = 85;
% % w_glc_enrich = 89;
% % %
% % w_idxs1 = w_glc_enrich;
% % w_idxs2 = [w_glc_enrich, w_glc_ec];
% % w_idxs3 = [w_glc_enrich, w_glt];
% % w_idxs4 = [w_glc_enrich, w_glc_ec, w_glt];
% % w_idxs5 = [w_glc_enrich, w_glc_ec, w_glt, ...
% %     w_treRates_glk_pgi,...
% %     w_g1p, w_udpg, w_t6p_tre, w_tre_ic, w_tre_ec];
% 
warray_cell = cell(ntests,1);
for i = 1:ntests % base,
    warray_cell{i} = blankWeight;
end
% 
warray_cell{1}(w_simHXK) = ones;
warray_cell{2}(w_simHXK) = ones;
warray_cell{3}(w_simHXK) = ones;
warray_cell{4}(w_simHXK) = ones;
warray_cell{5}(w_simHXK) = ones;
warray_cell{6}(w_simHXK) = ones;
warray_cell{7}(w_simHXK) = ones;
warray_cell{8}(w_simHXK) = ones;
% 
warray_cell{9}(w_simHXK) = ones;
warray_cell{10}(w_simHXK) = ones;
warray_cell{11}(w_simHXK) = ones;
warray_cell{12}(w_simHXK) = ones;
warray_cell{13}(w_simHXK) = ones;
warray_cell{14}(w_simHXK) = ones;
warray_cell{15}(w_simHXK) = ones;
warray_cell{16}(w_simHXK) = ones;
% 
% % previous
% warray_cell{1}(w_refY3M1) = ones * 1E-3;
% warray_cell{2}(w_refY3M1) = ones * 1E-2;
% warray_cell{3}(w_refY3M1) = ones * 1E-1; % 
% warray_cell{4}(w_refY3M1) = ones * 1E0; %
% warray_cell{5}(w_refY3M1) = ones * 3E0;
% warray_cell{6}(w_refY3M1) = ones * 1E1;
% warray_cell{7}(w_refY3M1) = ones * 1E2;
% warray_cell{8}(w_refY3M1) = ones * 1E3;
% 
% current
warray_cell{1}(w_refY3M1) = ones * 3E-3;
warray_cell{2}(w_refY3M1) = ones * 3E-2;
% logspace(-1,0,10)'
% 
warray_cell{3}(w_refY3M1) = 0.1000;
warray_cell{4}(w_refY3M1) = 0.1292;
warray_cell{5}(w_refY3M1) = 0.1668;
warray_cell{6}(w_refY3M1) = 0.2154;
warray_cell{7}(w_refY3M1) = 0.2783;
% 
warray_cell{8}(w_refY3M1) = 0.3594;
warray_cell{9}(w_refY3M1) = 0.4642;
warray_cell{10}(w_refY3M1) = 0.5995;
warray_cell{11}(w_refY3M1) = 0.7743;
warray_cell{12}(w_refY3M1) = 1.0000;
% 
warray_cell{13}(w_refY3M1) = ones * 3E0;
warray_cell{14}(w_refY3M1) = ones * 3E1;
warray_cell{15}(w_refY3M1) = ones * 3E2;
warray_cell{16}(w_refY3M1) = ones * 3E3;


% %% Parameter estimation
% % 
% selPars = pC16a;
% setup.selPars = selPars;
% % 
% cluster = parcluster('local');
% % pool = parpool(cluster,4);
% pool = parpool(cluster,16);
% parfor o = 1:16
% 
%     % disp
%     fprintf('o = %f.\n', o) %disp(o)
% 
%     % find in the precreated arrays
%     warray = warray_cell{o}; 
%     saveName = names_cell{o};
%     
%     % core options:
%     plength = length(selPars);
%     lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
%     ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;
%     x_temp = x(selPars);
% 
% %     % %% run check
% %     [error]=costfunSystemY3M1_FF_pRun(x_temp,canelas_SS,setup,xarray,data,dataset,NumberCycles,IC,selPars,warray);
% %     errorAnalysis_Y3M1;
% %     % %%
% 
% %     %
% %     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x_iter,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% %     plotMode = 2; % multiple simulations
% %     % clear selResCell
% %     selResCell{1}.T_FF01 = T_FF01_1;
% %     selResCell{1}.Y_FF01 = Y_FF01_1;
% %     selResCell{1}.V_FF01 = V_FF01_1;
% %     % 
% %     referencePlotSimulations_enrichment
% %     plotMode = 0;
% %     %
%         
%         % parameter estimation
% %         NumberCycles = 20;
%         sprintf('run num.%d tested',o)
%         tic
% %         [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun,x_temp,lb,ub,options,canelas_SS,setup,x_iter,data,dataset,NumberCycles,IC,selPars,warray);
%         [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun_enrichment,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,NumberCycles,IC,selPars,warray);
%         t = toc; 
%         disp(xres);
% 
%         % saving
%         parsave_Y3M2_cluster(saveName, xres, resnorm, residual, exitflag, warray, t, selPars);
% % % % %         save(saveName, 'xres', 'resnorm', 'residual', 'exitflag', 'warray', 't', 'selPars');
% end
% delete(pool)
% quit




%% DATA RECALL (not the multistart cases yet)
% xFullScale_99b_parCombs5_warray8.mat is lacking at start
% 
xAll1 = x;
%
namesHits1 = cell(1,1); namesHits1{1} = 'initial';
%
for i = 1:ntests
    % pC1_wA1
    temp_idx = i;
    loadName = sprintf(names_cell{temp_idx}, i);
    if exist(loadName,'file') == 2 % if exist
        load(loadName);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1 = [xAll1; x3];
        namesHits1 = [namesHits1; loadName];
    end
end
% 

% %% simulate results
% % 
% simRes1 = cell(1,length(namesHits1));
% for i = 1:length(namesHits1)
%     disp(i)
%     xSel = xAll1(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% % 
% 
% % %% select main case and save
%        
% % %%
% simRes1_16a2 = simRes1;
% save('tempRes_16a2.mat','simRes1_16a2')
%%
load('tempRes_16a2.mat','simRes1_16a2')


%% Visualizing results: simulation

% 
selResCell = simRes1_16a2;
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;

% % visualize (1:8), v_simHXK = ones, w_refY3M1 = increasing
% selResCell = simRes1_16a([1,2:9]);
% plotMode = 2;
% referencePlotSimulations_enrichment
% plotMode = 0;
% 
% % visualize (9:16), v_simHXK = increasing, w_refY3M1 = ones
% selResCell = simRes1_16a([1,10:17]);
% plotMode = 2;
% referencePlotSimulations_enrichment
% plotMode = 0;
% 


%% Visualizing results: parameters
nPoints = 16;
%
pE16a_w_simHXK_16a2 = ones(1,nPoints);
pE16a_w_refY3M1_16a2 = [3E-3,2E-2,...
                        0.1000,0.1292,0.1668,0.2154,0.2783,...
                        0.3594,0.4642,0.5995,0.7743,1.0000,...
                        3E0,3E1,3E2,3E3];
%


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


%% pE16a_18
% pE16a_18
% clf(1001)
fig1001 = figure(1001);
fig1001.Position = [1987 514 1024 420];

% semilogx 
subplot(1,2,1)
yyaxis left
semilogx(pE16a_w_refY3M1_16a2, pE16a2_error_simHXK, 'o-')
hold on
yyaxis right
semilogx(pE16a_w_refY3M1_16a2, pE16a2_error_pars, 'o-')
title('w_{HXK} = ones -> w_{refY3M1} (semilogx)')
xlim([min(pE16a_w_refY3M1_16a2) max(pE16a_w_refY3M1_16a2)])
% highlight point
% temp_id = 4;
temp_id = 10;
semilogx(pE16a_w_refY3M1_16a2(temp_id), pE16a2_error_simHXK(temp_id), 'o-', 'MarkerFaceColor','k')
semilogx(pE16a_w_refY3M1_16a2(temp_id), pE16a2_error_pars(temp_id), 'o-', 'MarkerFaceColor','k')
% 
legend('error simHXK','error parameters')
% 
hold off


% R roc-like curve
subplot(1,2,2)
plot(pE16a2_error_pars, pE16a2_error_simHXK, ...
    'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
xlabel('error_{pars}')
ylabel('error_{simHXK}')
title('roc-like curve')
hold off


%% Visualizing results: simulation. Checkign specific cases
% selResCell = simRes1_16a2([1,2:5]);
selResCell = simRes1_16a2([1,2:11]);
% selResCell = simRes1_16a2([1,2:7]);
% selResCell = simRes1_16a2([1,7:11]);
% selResCell = simRes1_16a2([1,2:10]);
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;


%% comparing changes in parameters using figure manuscript

% layout
maxWidth = 5;
maxRange = 1;
nSel = 4;
colRange = cool(nSel);
xSel = xAll1([1 4 11 17],:);
% xSel = xAll1([1 4 5 17],:);
% mBlue = [0, 0.4470, 0.7410];
% mOrange = [0.8500, 0.3250, 0.0980];

%
xlen = 1:length(parsGLK);
tempLabel = cell(length(parsGLK),1);
clf(103)
fh103 = figure(103);
for i = xlen
    wid = maxWidth * (abs(xSel(3,parsGLK(i)) - xSel(4,parsGLK(i)))/maxRange);
    line([xlen(i) xlen(i)],[min([xSel(3,parsGLK(i)) xSel(4,parsGLK(i))]) max([xSel(3,parsGLK(i)) xSel(4,parsGLK(i))])],...
        'Color','black','LineWidth',wid)
    tempLabel{i} = strtok(legenda.parameters{parsGLK(i)}(9:end),', ');
end
hold on
%
for i = 1:nSel
    plot(xlen,abs(xSel(i,parsGLK)),'color',colRange(i,:),'LineStyle',':')
    plot(xlen,-abs(xSel(i,parsGLK)),'color',colRange(i,:),'LineStyle',':')
end
%
h = zeros(1,nSel);
for i = 1:nSel
    h(i) = scatter(xlen,xSel(i,parsGLK),'filled', ...
        'MarkerEdgeColor', colRange(i,:), ...
        'MarkerFaceColor', colRange(i,:));
end
xticks(xlen)
xticklabels(tempLabel)
ylim([-3 3])
legend(h,'xY3M2','x4','x11','xY3M1',...
    'Location','SouthOutside','Orientation','horizontal')
title('parsGLK')
hold off


%%
% % previous
% xFinal = xAll1(15,:);
% xFinal(30) = xAll1(14,30); 
% xTemp = [xAll1([1 15],:); xFinal];
% 16a2
xFinal = xAll1(11,:);
% xFinal(30) = xAll1(14,30); 
xTemp = [xAll1([1 17],:); xFinal];
%%
simRes1_temp = cell(1,3);
for i = 1:3
    disp(i)
    xSel = xTemp(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_temp{i}.T_FF01 = T_FF01;
    simRes1_temp{i}.Y_FF01 = Y_FF01;
    simRes1_temp{i}.V_FF01 = V_FF01;
end
selResCell = simRes1_temp;
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;
%%
x16a2 = xFinal;
save('x16a2.mat','x16a2','xFinal')


%% parameter comparison
% 
load('x16a.mat','x16a')
%% 
xFinal_16a1 = x16a;
xFinal_16a2 = x16a2;
%%
temp_pars = [32    33    34];
xtemp2 =[xAll1(17,temp_pars);...
    xFinal_16a1(temp_pars);...
    xFinal_16a2(temp_pars)]
%%
% p.HXK1_Kglc=10.^x(32).*0.08; %(UNIT!)
% p.HXK1_Kt6p=10.^x(33).*0.2; %(UNIT!)
% p.HXK1_kcat=.2*10.^x(34).*4.75; %(UNIT!) 

temp_pars_real = [0.08 .* 10 .^ xtemp2(:,1),...
    0.2 .* 10 .^ xtemp2(:,2),...
    0.2 .* 4.75 .* 10 .^ xtemp2(:,3)]

% fold change
temp_pars_fold = [ temp_pars_real(2,:) ./ temp_pars_real(1,:); ...
                    temp_pars_real(3,:) ./ temp_pars_real(1,:)]



















