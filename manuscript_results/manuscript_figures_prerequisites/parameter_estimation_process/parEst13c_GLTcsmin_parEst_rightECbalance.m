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


%% No enrichment simulations
% % % % % plotflag = 2; % variables by iteration
% % % % choosedataset
% % % % 
% % % % setup.experiment = 1;
% % % % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % plotMode = 2; % multiple simulations
% % % % % clear selResCell
% % % % selResCell{1}.T_FF01 = T_FF01_1;
% % % % selResCell{1}.Y_FF01 = Y_FF01_1;
% % % % selResCell{1}.V_FF01 = V_FF01_1;
% % % % 
% % % % referencePlotSimulations
% % % % plotMode = 0;
% % % % 
% % % % setup.clamp_GLCec = 0;
% % % % % toc


% % %% Simulating enrichment
% % legendaMetabolites_addEnrichment;
% % % plotflag = 2; % variables by iteration
% % choosedataset
% % setup.experiment = 1;
% % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % plotMode = 2; % multiple simulations
% % % clear selResCell
% % selResCell{1}.T_FF01 = T_FF01_1;
% % selResCell{1}.Y_FF01 = Y_FF01_1;
% % selResCell{1}.V_FF01 = V_FF01_1;
% % % 
% % referencePlotSimulations_enrichment
% % plotMode = 0;
% % 
% % % setup.clamp_GLCec = 0;
% % % toc


%% PSA on this initial simulations
x2 = x; x2(164) = 5;
plotMode = 0;
choosedataset
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
[T_FF01_2,Y_FF01_2,V_FF01_2] = simulate_FF_enrichment(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
% clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;
selResCell{2}.T_FF01 = T_FF01_2;
selResCell{2}.Y_FF01 = Y_FF01_2;
selResCell{2}.V_FF01 = V_FF01_2;
% 
referencePlotSimulations_enrichment




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Parameter estimation setup

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
ntests = 15;
nMS = 20;

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
names_cell{1} = 'FF_pE13c_pC1_wA1_nS%d.mat';
names_cell{2} = 'FF_pE13c_pC1_wA2_nS%d.mat';
names_cell{3} = 'FF_pE13c_pC1_wA3_nS%d.mat';
names_cell{4} = 'FF_pE13c_pC1_wA4_nS%d.mat';
names_cell{5} = 'FF_pE13c_pC1_wA5_nS%d.mat';
% 
names_cell{6} = 'FF_pE13c_pC2_wA1_nS%d.mat';
names_cell{7} = 'FF_pE13c_pC2_wA2_nS%d.mat';
names_cell{8} = 'FF_pE13c_pC2_wA3_nS%d.mat';
names_cell{9} = 'FF_pE13c_pC2_wA4_nS%d.mat';
names_cell{10} = 'FF_pE13c_pC2_wA5_nS%d.mat';
% 
names_cell{11} = 'FF_pE13c_pC3_wA1_nS%d.mat';
names_cell{12} = 'FF_pE13c_pC3_wA2_nS%d.mat';
names_cell{13} = 'FF_pE13c_pC3_wA3_nS%d.mat';
names_cell{14} = 'FF_pE13c_pC3_wA4_nS%d.mat';
names_cell{15} = 'FF_pE13c_pC3_wA5_nS%d.mat';
%
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

pC1 = [parsGLT, parsATH1ec];
pC2 = [parsGLT, parsHXK, parsTreExport];
pC3 = [parsGLT, parsHXK, parsTreCycle, parsTreExport];
%
parComb_cell = cell(ntests,1);
% 
parComb_cell{1} = pC1;
parComb_cell{2} = pC1;
parComb_cell{3} = pC1;
parComb_cell{4} = pC1;
parComb_cell{5} = pC1;
% 
parComb_cell{6} = pC2;
parComb_cell{7} = pC2;
parComb_cell{8} = pC2;
parComb_cell{9} = pC2;
parComb_cell{10} = pC2;
% 
parComb_cell{11} = pC3;
parComb_cell{12} = pC3;
parComb_cell{13} = pC3;
parComb_cell{14} = pC3;
parComb_cell{15} = pC3;
% 

% weights:
w_glc_ic = 6;
w_glc_ec = 36;
w_g1p = 21;
w_udpg = 24;
w_t6p_tre = [26 25];
w_treRates_glk_pgi = [55 59 57 58 56 40 41];
w_glt = 39;
w_gly = [86 87 88];
w_tre_ic = 25;
w_tre_ec = 37;
w_ath1 = [83 84];
    w_ath1_ec = 83;
    w_ath1_vac = 84; 
w_agt1 = 85;
w_glc_enrich = 89;
%
w_idxs1 = w_glc_enrich;
w_idxs2 = [w_glc_enrich, w_glc_ec];
w_idxs3 = [w_glc_enrich, w_glt];
w_idxs4 = [w_glc_enrich, w_glc_ec, w_glt];
w_idxs5 = [w_glc_enrich, w_glc_ec, w_glt, ...
    w_treRates_glk_pgi,...
    w_g1p, w_udpg, w_t6p_tre, w_tre_ic, w_tre_ec];
%
warray_cell = cell(ntests,1);
for i = 1:ntests % base,
    warray_cell{i} = blankWeight;
end
% 
warray_cell{1}(w_idxs1) = ones;
warray_cell{2}(w_idxs2) = ones;
warray_cell{3}(w_idxs3) = ones;
warray_cell{4}(w_idxs4) = ones;
warray_cell{5}(w_idxs5) = ones;
% 
warray_cell{6}(w_idxs1) = ones;
warray_cell{7}(w_idxs2) = ones;
warray_cell{8}(w_idxs3) = ones;
warray_cell{9}(w_idxs4) = ones;
warray_cell{10}(w_idxs5) = ones;
% 
warray_cell{11}(w_idxs1) = ones;
warray_cell{12}(w_idxs2) = ones;
warray_cell{13}(w_idxs3) = ones;
warray_cell{14}(w_idxs4) = ones;
warray_cell{15}(w_idxs5) = ones;
% 


%% Selected runs
% 1*165
selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS
selectedRuns_idx = selCase
% if selCase <= 24
%     selectedRuns_idx = [1,17,33] + 1 * (selCase - 1)
% else
%     selectedRuns_idx = [1:15] + 15 * (selCase - 25) + 384;
% end
% disp(selectedRuns_idx)
% % % % selectedRuns_idx = 1:ntests;
% 
% 
%% PARAMETER ESTIMATION: NO MS
for o = selectedRuns_idx
    
    for u = 1:nMS
        % disp
        fprintf('o = %f, u = %f.\n', o, u) %disp(o)

        % find in the precreated arrays
        selPars = parComb_cell{o};
        warray = warray_cell{o}; 
        tempName = names_cell{o};
        saveName = sprintf(tempName,u);

        % core options:
        x_iter = x;
        if u == 1 % initial
        elseif u == 2 % ath1ec hydrolysis
            x_iter(164) = 0;
        end
        x_temp = x_iter(selPars);
        if u > 2 % rand
            rng(u); x_temp = x_temp - 1 + 2 * rand(size(x_temp));%x_iter
        end
        plength = length(selPars);
        lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
        ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;
        if(o >= 6)
            lb(end-12:end) = -20;
            ub(end-12:end) = 20;
        else
            lb(end-2:end) = -20;
            ub(end-2:end) = 20;
        end

    %     % %% run check
    %     [error]=costfunSystemY3M1_FF_pRun(x_temp,canelas_SS,setup,xarray,data,dataset,NumberCycles,IC,selPars,warray);
    %     errorAnalysis_Y3M1;
    %     % %%

%         %
%         [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x_iter,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%         plotMode = 2; % multiple simulations
%         % clear selResCell
%         selResCell{1}.T_FF01 = T_FF01_1;
%         selResCell{1}.Y_FF01 = Y_FF01_1;
%         selResCell{1}.V_FF01 = V_FF01_1;
%         % 
%         referencePlotSimulations_enrichment
%         plotMode = 0;
%         %
        
        % parameter estimation
        sprintf('run num.%d tested',o)
        tic
%         [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun,x_temp,lb,ub,options,canelas_SS,setup,x_iter,data,dataset,NumberCycles,IC,selPars,warray);
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun_enrichment,x_temp,lb,ub,options,canelas_SS,setup,x_iter,data,dataset,NumberCycles,IC,selPars,warray);
        t = toc; 
        disp(xres);

        % saving
    % % % %     parsave_Y3M2_cluster(saveName, xres, resnorm, residual, exitflag, warray, t, selPars);
        save(saveName, 'xres', 'resnorm', 'residual', 'exitflag', 'warray', 't', 'selPars');
    end
    
end
% disp('parallel run below')
% % save
% saveName = 'FF_pE4_x1.mat';
% save(saveName,'allRes_xres','allRes_resnorm','allRes_residual','allRes_exitflag','allRes_t','allRes_warray','allRes_selPars');
% % % % delete(pool)
% quit


%% DATA RECALL (not the multistart cases yet)
% xFullScale_99b_parCombs5_warray8.mat is lacking at start
% 
xAll1_pC1_wA1 = x;
xAll1_pC1_wA2 = x;
xAll1_pC1_wA3 = x;
xAll1_pC1_wA4 = x;
xAll1_pC1_wA5 = x;
% 
xAll1_pC2_wA1 = x;
xAll1_pC2_wA2 = x;
xAll1_pC2_wA3 = x;
xAll1_pC2_wA4 = x;
xAll1_pC2_wA5 = x;
% 
xAll1_pC3_wA1 = x;
xAll1_pC3_wA2 = x;
xAll1_pC3_wA3 = x;
xAll1_pC3_wA4 = x;
xAll1_pC3_wA5 = x;
%
namesHits1_pC1_wA1 = cell(1,1); namesHits1_pC1_wA1{1} = 'initial';
namesHits1_pC1_wA2 = cell(1,1); namesHits1_pC1_wA2{1} = 'initial';
namesHits1_pC1_wA3 = cell(1,1); namesHits1_pC1_wA3{1} = 'initial';
namesHits1_pC1_wA4 = cell(1,1); namesHits1_pC1_wA4{1} = 'initial';
namesHits1_pC1_wA5 = cell(1,1); namesHits1_pC1_wA5{1} = 'initial';
%
namesHits1_pC2_wA1 = cell(1,1); namesHits1_pC2_wA1{1} = 'initial';
namesHits1_pC2_wA2 = cell(1,1); namesHits1_pC2_wA2{1} = 'initial';
namesHits1_pC2_wA3 = cell(1,1); namesHits1_pC2_wA3{1} = 'initial';
namesHits1_pC2_wA4 = cell(1,1); namesHits1_pC2_wA4{1} = 'initial';
namesHits1_pC2_wA5 = cell(1,1); namesHits1_pC2_wA5{1} = 'initial';
%
namesHits1_pC3_wA1 = cell(1,1); namesHits1_pC3_wA1{1} = 'initial';
namesHits1_pC3_wA2 = cell(1,1); namesHits1_pC3_wA2{1} = 'initial';
namesHits1_pC3_wA3 = cell(1,1); namesHits1_pC3_wA3{1} = 'initial';
namesHits1_pC3_wA4 = cell(1,1); namesHits1_pC3_wA4{1} = 'initial';
namesHits1_pC3_wA5 = cell(1,1); namesHits1_pC3_wA5{1} = 'initial';
%
for i = 1:nMS
    % pC1_wA1
    temp_idx = 1;
    loadName_pC1_wA1 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC1_wA1,'file') == 2 % if exist
        load(loadName_pC1_wA1);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC1_wA1 = [xAll1_pC1_wA1; x3];
        namesHits1_pC1_wA1 = [namesHits1_pC1_wA1; loadName_pC1_wA1];
    end
    % pC1_wA2
    temp_idx = 2;
    loadName_pC1_wA2 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC1_wA2,'file') == 2 % if exist
        load(loadName_pC1_wA2);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC1_wA2 = [xAll1_pC1_wA2; x3];
        namesHits1_pC1_wA2 = [namesHits1_pC1_wA2; loadName_pC1_wA2];
    end
    % pC1_wA3
    temp_idx = 3;
    loadName_pC1_wA3 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC1_wA3,'file') == 2 % if exist
        load(loadName_pC1_wA3);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC1_wA3 = [xAll1_pC1_wA3; x3];
        namesHits1_pC1_wA3 = [namesHits1_pC1_wA3; loadName_pC1_wA3];
    end
    % pC1_wA4
    temp_idx = 4;
    loadName_pC1_wA4 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC1_wA4,'file') == 2 % if exist
        load(loadName_pC1_wA4);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC1_wA4 = [xAll1_pC1_wA4; x3];
        namesHits1_pC1_wA4 = [namesHits1_pC1_wA4; loadName_pC1_wA4];
    end
    % pC1_wA5
    temp_idx = 5;
    loadName_pC1_wA5 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC1_wA5,'file') == 2 % if exist
        load(loadName_pC1_wA5);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC1_wA5 = [xAll1_pC1_wA5; x3];
        namesHits1_pC1_wA5 = [namesHits1_pC1_wA5; loadName_pC1_wA5];
    end
end
%
for i = 1:nMS
    % pC2_wA1
    temp_idx = 6;
    loadName_pC2_wA1 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC2_wA1,'file') == 2 % if exist
        load(loadName_pC2_wA1);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC2_wA1 = [xAll1_pC2_wA1; x3];
        namesHits1_pC2_wA1 = [namesHits1_pC2_wA1; loadName_pC2_wA1];
    end
    % pC2_wA2
    temp_idx = 7;
    loadName_pC2_wA2 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC2_wA2,'file') == 2 % if exist
        load(loadName_pC2_wA2);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC2_wA2 = [xAll1_pC2_wA2; x3];
        namesHits1_pC2_wA2 = [namesHits1_pC2_wA2; loadName_pC2_wA2];
    end
    % pC2_wA3
    temp_idx = 8;
    loadName_pC2_wA3 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC2_wA3,'file') == 2 % if exist
        load(loadName_pC2_wA3);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC2_wA3 = [xAll1_pC2_wA3; x3];
        namesHits1_pC2_wA3 = [namesHits1_pC2_wA3; loadName_pC2_wA3];
    end
    % pC2_wA4
    temp_idx = 9;
    loadName_pC2_wA4 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC2_wA4,'file') == 2 % if exist
        load(loadName_pC2_wA4);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC2_wA4 = [xAll1_pC2_wA4; x3];
        namesHits1_pC2_wA4 = [namesHits1_pC2_wA4; loadName_pC2_wA4];
    end
    % pC2_wA5
    temp_idx = 10;
    loadName_pC2_wA5 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC2_wA5,'file') == 2 % if exist
        load(loadName_pC2_wA5);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC2_wA5 = [xAll1_pC2_wA5; x3];
        namesHits1_pC2_wA5 = [namesHits1_pC2_wA5; loadName_pC2_wA5];
    end
end
%
for i = 1:nMS
    % pC3_wA1
    temp_idx = 11;
    loadName_pC3_wA1 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC3_wA1,'file') == 2 % if exist
        load(loadName_pC3_wA1);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC3_wA1 = [xAll1_pC3_wA1; x3];
        namesHits1_pC3_wA1 = [namesHits1_pC3_wA1; loadName_pC3_wA1];
    end
    % pC3_wA2
    temp_idx = 12;
    loadName_pC3_wA2 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC3_wA2,'file') == 2 % if exist
        load(loadName_pC3_wA2);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC3_wA2 = [xAll1_pC3_wA2; x3];
        namesHits1_pC3_wA2 = [namesHits1_pC3_wA2; loadName_pC3_wA2];
    end
    % pC3_wA3
    temp_idx = 13;
    loadName_pC3_wA3 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC3_wA3,'file') == 2 % if exist
        load(loadName_pC3_wA3);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC3_wA3 = [xAll1_pC3_wA3; x3];
        namesHits1_pC3_wA3 = [namesHits1_pC3_wA3; loadName_pC3_wA3];
    end
    % pC3_wA4
    temp_idx = 14;
    loadName_pC3_wA4 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC3_wA4,'file') == 2 % if exist
        load(loadName_pC3_wA4);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC3_wA4 = [xAll1_pC3_wA4; x3];
        namesHits1_pC3_wA4 = [namesHits1_pC3_wA4; loadName_pC3_wA4];
    end
    % pC3_wA5
    temp_idx = 15;
    loadName_pC3_wA5 = sprintf(names_cell{temp_idx}, i);
    if exist(loadName_pC3_wA5,'file') == 2 % if exist
        load(loadName_pC3_wA5);
        selPars = parComb_cell{temp_idx};
        x3 = x;
        x3(selPars) = xres;
        xAll1_pC3_wA5 = [xAll1_pC3_wA5; x3];
        namesHits1_pC3_wA5 = [namesHits1_pC3_wA5; loadName_pC3_wA5];
    end
end


%% simulate results
% 
simRes1_pC1_wA1 = cell(1,length(namesHits1_pC1_wA1));
simRes1_pC1_wA2 = cell(1,length(namesHits1_pC1_wA2));
simRes1_pC1_wA3 = cell(1,length(namesHits1_pC1_wA3));
simRes1_pC1_wA4 = cell(1,length(namesHits1_pC1_wA4));
simRes1_pC1_wA5 = cell(1,length(namesHits1_pC1_wA5));
% 
simRes1_pC2_wA1 = cell(1,length(namesHits1_pC2_wA1));
simRes1_pC2_wA2 = cell(1,length(namesHits1_pC2_wA2));
simRes1_pC2_wA3 = cell(1,length(namesHits1_pC2_wA3));
simRes1_pC2_wA4 = cell(1,length(namesHits1_pC2_wA4));
simRes1_pC2_wA5 = cell(1,length(namesHits1_pC2_wA5));
% 
simRes1_pC3_wA1 = cell(1,length(namesHits1_pC3_wA1));
simRes1_pC3_wA2 = cell(1,length(namesHits1_pC3_wA2));
simRes1_pC3_wA3 = cell(1,length(namesHits1_pC3_wA3));
simRes1_pC3_wA4 = cell(1,length(namesHits1_pC3_wA4));
simRes1_pC3_wA5 = cell(1,length(namesHits1_pC3_wA5));


cluster = parcluster('local');

% pC1
pool = parpool(cluster,4);
parfor i = 1:length(namesHits1_pC1_wA1)
    disp(i)
    xSel = xAll1_pC1_wA1(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC1_wA1{i}.T_FF01 = T_FF01;
    simRes1_pC1_wA1{i}.Y_FF01 = Y_FF01;
    simRes1_pC1_wA1{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC1_wA2)
    disp(i)
    xSel = xAll1_pC1_wA2(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC1_wA2{i}.T_FF01 = T_FF01;
    simRes1_pC1_wA2{i}.Y_FF01 = Y_FF01;
    simRes1_pC1_wA2{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC1_wA3)
    disp(i)
    xSel = xAll1_pC1_wA3(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC1_wA3{i}.T_FF01 = T_FF01;
    simRes1_pC1_wA3{i}.Y_FF01 = Y_FF01;
    simRes1_pC1_wA3{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC1_wA4)
    disp(i)
    xSel = xAll1_pC1_wA4(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC1_wA4{i}.T_FF01 = T_FF01;
    simRes1_pC1_wA4{i}.Y_FF01 = Y_FF01;
    simRes1_pC1_wA4{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC1_wA5)
    disp(i)
    xSel = xAll1_pC1_wA5(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC1_wA5{i}.T_FF01 = T_FF01;
    simRes1_pC1_wA5{i}.Y_FF01 = Y_FF01;
    simRes1_pC1_wA5{i}.V_FF01 = V_FF01;
end
%

% pC2
parfor i = 1:length(namesHits1_pC2_wA1)
    disp(i)
    xSel = xAll1_pC2_wA1(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC2_wA1{i}.T_FF01 = T_FF01;
    simRes1_pC2_wA1{i}.Y_FF01 = Y_FF01;
    simRes1_pC2_wA1{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC2_wA2)
    disp(i)
    xSel = xAll1_pC2_wA2(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC2_wA2{i}.T_FF01 = T_FF01;
    simRes1_pC2_wA2{i}.Y_FF01 = Y_FF01;
    simRes1_pC2_wA2{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC2_wA3)
    disp(i)
    xSel = xAll1_pC2_wA3(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC2_wA3{i}.T_FF01 = T_FF01;
    simRes1_pC2_wA3{i}.Y_FF01 = Y_FF01;
    simRes1_pC2_wA3{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC2_wA4)
    disp(i)
    xSel = xAll1_pC2_wA4(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC2_wA4{i}.T_FF01 = T_FF01;
    simRes1_pC2_wA4{i}.Y_FF01 = Y_FF01;
    simRes1_pC2_wA4{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC2_wA5)
    disp(i)
    xSel = xAll1_pC2_wA5(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC2_wA5{i}.T_FF01 = T_FF01;
    simRes1_pC2_wA5{i}.Y_FF01 = Y_FF01;
    simRes1_pC2_wA5{i}.V_FF01 = V_FF01;
end
%

% pC3
parfor i = 1:length(namesHits1_pC3_wA1)
    disp(i)
    xSel = xAll1_pC3_wA1(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC3_wA1{i}.T_FF01 = T_FF01;
    simRes1_pC3_wA1{i}.Y_FF01 = Y_FF01;
    simRes1_pC3_wA1{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC3_wA2)
    disp(i)
    xSel = xAll1_pC3_wA2(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC3_wA2{i}.T_FF01 = T_FF01;
    simRes1_pC3_wA2{i}.Y_FF01 = Y_FF01;
    simRes1_pC3_wA2{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC3_wA3)
    disp(i)
    xSel = xAll1_pC3_wA3(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC3_wA3{i}.T_FF01 = T_FF01;
    simRes1_pC3_wA3{i}.Y_FF01 = Y_FF01;
    simRes1_pC3_wA3{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC3_wA4)
    disp(i)
    xSel = xAll1_pC3_wA4(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC3_wA4{i}.T_FF01 = T_FF01;
    simRes1_pC3_wA4{i}.Y_FF01 = Y_FF01;
    simRes1_pC3_wA4{i}.V_FF01 = V_FF01;
end
%
parfor i = 1:length(namesHits1_pC3_wA5)
    disp(i)
    xSel = xAll1_pC3_wA5(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1_pC3_wA5{i}.T_FF01 = T_FF01;
    simRes1_pC3_wA5{i}.Y_FF01 = Y_FF01;
    simRes1_pC3_wA5{i}.V_FF01 = V_FF01;
end
%

delete(pool)

%% select main case and save
% %
simRes1_pC1 = [simRes1_pC1_wA1(1:2),...
                simRes1_pC1_wA2(2),...
                simRes1_pC1_wA3(2),...
                simRes1_pC1_wA4(2),...
                simRes1_pC1_wA5(2)];
%
namesHits1_pC1 = [namesHits1_pC1_wA1(1:2)', ...
                namesHits1_pC1_wA2(2), ...
                namesHits1_pC1_wA3(2), ...
                namesHits1_pC1_wA4(2), ...
                namesHits1_pC1_wA5(2)];

% %
simRes1_pC2 = [simRes1_pC2_wA1(1:2),...
                simRes1_pC2_wA2(2),...
                simRes1_pC2_wA3(2),...
                simRes1_pC2_wA4(2),...
                simRes1_pC2_wA5(2)];
%
namesHits1_pC2 = [namesHits1_pC2_wA1(1:2)', ...
                namesHits1_pC2_wA2(2), ...
                namesHits1_pC2_wA3(2), ...
                namesHits1_pC2_wA4(2), ...
                namesHits1_pC2_wA5(2)];
            
% %
simRes1_pC3 = [simRes1_pC3_wA1(1:2),...
                simRes1_pC3_wA2(2),...
                simRes1_pC3_wA3(2),...
                simRes1_pC3_wA4(2),...
                simRes1_pC3_wA5(2)];
%
namesHits1_pC3 = [namesHits1_pC3_wA1(1:2)', ...
                namesHits1_pC3_wA2(2), ...
                namesHits1_pC3_wA3(2), ...
                namesHits1_pC3_wA4(2), ...
                namesHits1_pC3_wA5(2)];
            
% %
simRes1_main = [simRes1_pC1,...
                simRes1_pC2(2:end),...
                simRes1_pC3(2:end)];
%
namesHits1_main = [namesHits1_pC1, ...
                namesHits1_pC2(2:end), ...
                namesHits1_pC3(2:end)];
            
            
%%
save('tempRes_13c.mat',...
    'simRes1_pC1_wA1', 'simRes1_pC1_wA2', 'simRes1_pC1_wA3', 'simRes1_pC1_wA4', 'simRes1_pC1_wA5', ...
    'simRes1_pC2_wA1', 'simRes1_pC2_wA2', 'simRes1_pC2_wA3', 'simRes1_pC2_wA4', 'simRes1_pC2_wA5', ...
    'simRes1_pC3_wA1', 'simRes1_pC3_wA2', 'simRes1_pC3_wA3', 'simRes1_pC3_wA4', 'simRes1_pC3_wA5', ...
    'simRes1_pC1', 'simRes1_pC2', 'simRes1_pC3', 'simRes1_main')


%%
load('tempRes_13c.mat',...
    'simRes1_pC1_wA1', 'simRes1_pC1_wA2', 'simRes1_pC1_wA3', 'simRes1_pC1_wA4', 'simRes1_pC1_wA5', ...
    'simRes1_pC2_wA1', 'simRes1_pC2_wA2', 'simRes1_pC2_wA3', 'simRes1_pC2_wA4', 'simRes1_pC2_wA5', ...
    'simRes1_pC3_wA1', 'simRes1_pC3_wA2', 'simRes1_pC3_wA3', 'simRes1_pC3_wA4', 'simRes1_pC3_wA5', ...
    'simRes1_pC1', 'simRes1_pC2', 'simRes1_pC3', 'simRes1_main')


%% Visualizing results
% 
% select case to plot
selResCell = simRes1_pC1_wA1; % no enrichment
% selResCell = simRes1_pC1_wA2; % no enrichment
% selResCell = simRes1_pC1_wA3; % no enrichment. Better rates.
% selResCell = simRes1_pC1_wA4; % no enrichment. Better rates.
% selResCell = simRes1_pC1_wA5; % no enrichment. Better rates.
% 
% selResCell = simRes1_pC2_wA1; % no enrichment
% selResCell = simRes1_pC2_wA2; % no enrichment
% selResCell = simRes1_pC2_wA3; % no change
% selResCell = simRes1_pC2_wA4; % no change
% selResCell = simRes1_pC2_wA5; % no change
% 
% selResCell = simRes1_pC3_wA1;
% selResCell = simRes1_pC3_wA2;
% selResCell = simRes1_pC3_wA3;
% selResCell = simRes1_pC3_wA4;
% selResCell = simRes1_pC3_wA5;
% 
% selResCell = simRes1_pC1;
% selResCell = simRes1_pC2;
% selResCell = simRes1_pC3;
% selResCell = simRes1_main;

plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;


% % % % %% screening for correct enrichment level
% % % % simRes1_all = [simRes1_pC1_wA1, simRes1_pC1_wA2(2:end), simRes1_pC1_wA3(2:end), simRes1_pC1_wA4(2:end), simRes1_pC1_wA5(2:end), ...
% % % %                 simRes1_pC2_wA1(2:end), simRes1_pC2_wA2(2:end), simRes1_pC2_wA3(2:end), simRes1_pC2_wA4(2:end), simRes1_pC2_wA5(2:end), ...
% % % %                 simRes1_pC3_wA1(2:end), simRes1_pC3_wA2(2:end), simRes1_pC3_wA3(2:end), simRes1_pC3_wA4(2:end), simRes1_pC3_wA5(2:end)];
% % % % %%
% % % % n_temp = length(simRes1_all);
% % % % t_intend = 250;
% % % % idxs = [];
% % % % % figure
% % % % for i = 1:n_temp
% % % %     [~,idx] = min(abs(t_intend - simRes1_all{i}.T_FF01));
% % % %     enrich_temp = simRes1_all{i}.Y_FF01(idx,65)./simRes1_all{i}.Y_FF01(idx,36);
% % % %     % plot option
% % % % %     scatter(i, enrich_temp, 'filled','MarkerEdgeColor','black',...
% % % % %               'MarkerFaceColor','black'), hold on
% % % %     % select points
% % % %     if((enrich_temp > 0.5) && (enrich_temp < 0.63))
% % % %         idxs = [idxs, i];
% % % %     end
% % % % end
% % % % 
% % % % %% plot rather specific cases
% % % % % 
% % % % selResCell = simRes1_all([1, idxs]);
% % % % % 
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;
% % % % 
% % % % 
% % % % %% Screening for the specific case that looked better
% % % % n_temp = length(simRes1_pC1_wA2);
% % % % t_intend = 250;
% % % % idxs2 = [];
% % % % % figure
% % % % for i = 1:n_temp
% % % %     [~,idx] = min(abs(t_intend - simRes1_pC1_wA2{i}.T_FF01));
% % % %     enrich_temp = simRes1_pC1_wA2{i}.Y_FF01(idx,65)./simRes1_pC1_wA2{i}.Y_FF01(idx,36);
% % % %     % plot option
% % % % %     scatter(i, enrich_temp, 'filled','MarkerEdgeColor','black',...
% % % % %               'MarkerFaceColor','black'), hold on
% % % %     % select points
% % % %     if((enrich_temp > 0.5) && (enrich_temp < 0.63))
% % % %         idxs2 = [idxs2, i];
% % % %     end
% % % % end
% % % % 
% % % % 
% % % % %% plot rather specific cases
% % % % % 
% % % % selResCell = simRes1_pC1_wA2([1, idxs2]);
% % % % % 
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;
% % % % 
% % % % 
% % % % %% more in-depth study of the selected set
% % % % x_res = xAll1_pC1_wA2([1, idxs2],:);
% % % % % pC1 = [35    36    38   163   164   165];
% % % % 
% % % % %  
% % % % figure,
% % % % len_temp = length(pC1);
% % % % for i = 1:len_temp
% % % %     sp = subplot(2,3,i);
% % % %     histogram(x_res(2:end,pC1(i)))
% % % %     hold on
% % % %     line([x_res(1,pC1(i)) x_res(1,pC1(i))],sp.YLim,'Color','red','LineStyle','--')
% % % % end
% % % % 
% % % % %% replot unique 
% % % % selResCell = simRes1_pC1_wA2([1, 2]);
% % % % % 
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;
% % % % 
% % % % 
% % % % %% 
% % % % % 
% % % % xtest = x_res([1 1 1 2],:);
% % % % xtest(2,[35 36 38]) = xtest(4,[35 36 38]);
% % % % xtest(3,[163 164 165]) = xtest(4,[163 164 165]);
% % % % 
% % % % simRes1_test = cell(1,4);
% % % % for i = 1:4
% % % %     disp(i)
% % % %     xSel = xtest(i,:);
% % % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes1_test{i}.T_FF01 = T_FF01;
% % % %     simRes1_test{i}.Y_FF01 = Y_FF01;
% % % %     simRes1_test{i}.V_FF01 = V_FF01;
% % % % end
% % % % %%
% % % % selResCell = simRes1_test;
% % % % 
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;
% % % % 
% % % % 
% % % % %% selected cases to save
% % % % x13b_start = xtest(1,:);
% % % % x13b_GLTonly = xtest(2,:);
% % % % x13b_ATH1only = xtest(3,:);
% % % % x13b_end = xtest(4,:);
% % % % save('x13b.mat','x13b_start','x13b_GLTonly','x13b_ATH1only','x13b_end');
% % % % 
% % % % 
% % % % %%
% % % % % % 
% % % % % simRes1_test = cell(1,4);
% % % % for i = 4%[1,4,3]
% % % %     xSel = xtest(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % end
% % % % 
% % % % 
% % % % %% issue found. We'll re-optimize the setup again
% % % % setup.TRE_recirculation_rightBalances = 1;
% % % % 
% % % % xtest = x_res([1 2],:);
% % % % 
% % % % simRes1_test = cell(1,2);
% % % % for i = 1:2
% % % %     disp(i)
% % % %     xSel = xtest(i,:);
% % % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes1_test{i}.T_FF01 = T_FF01;
% % % %     simRes1_test{i}.Y_FF01 = Y_FF01;
% % % %     simRes1_test{i}.V_FF01 = V_FF01;
% % % % end
% % % % 
% % % % selResCell = simRes1_test;
% % % % 
% % % % plotMode = 2;
% % % % referencePlotSimulations_enrichment
% % % % plotMode = 0;
% % % % 
% % % % setup.TRE_recirculation_rightBalances = 0;





% %% testing area
% n_temp = length(simRes1_pC1_wA1);
% t_intend = 250;
% idxs2 = [];
% % figure
% for i = 1:n_temp
%     [~,idx] = min(abs(t_intend - simRes1_pC1_wA1{i}.T_FF01));
%     enrich_temp = simRes1_pC1_wA1{i}.Y_FF01(idx,65)./simRes1_pC1_wA1{i}.Y_FF01(idx,36);
%     % plot option
% %     scatter(i, enrich_temp, 'filled','MarkerEdgeColor','black',...
% %               'MarkerFaceColor','black'), hold on
%     % select points
%     if((enrich_temp > 0.5) && (enrich_temp < 0.58))
%         idxs2 = [idxs2, i];
%     end
% end
% 
% % %% plot rather specific cases
% % 
% selResCell = simRes1_pC1_wA1([1, idxs2]);
% % 
% plotMode = 2;
% referencePlotSimulations_enrichment
% plotMode = 0;



%% Understanding this specific case that worked better
% Which parameters change and in which direction?
% How does this alter the recirculation fluxes?




% % %%
% for i = 1:15
%     fid = fopen('parEst13c_GLTcsmin_parEst_rightECbalance.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d; % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE13c_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end

