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

%% (2021 09 17) Adjustment glk
% load('x16a.mat','x16a')
% x = x16a;
% x16c_E_start = x_test(1,:);
% x16c_E_final = x_test(2,:);
load('x16c_E_TPS2.mat', 'x16c_E_start', 'x16c_E_final')
x = x16c_E_final;

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
% % % 
% plotMode = 2;
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
warray_list = [1E-3 2E-3 5E-3,...
                1E-2 2E-2 5E-2,...
                1E-1 2E-1 5E-1,...
                1E0 2E0 5E0,...
                1E1 2E1 5E1,...
                1E2 2E2 5E2,...
                1E3 2E3 5E3];
ntests = length(warray_list);

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
for i = 1:ntests
    names_cell{i} = sprintf('FF_pE16d_TRE_wA%d.mat',i);
end
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
    parsTreCycle2 = [parsPGM1, parsTPS1, parsTPS2, parsNTH1];
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
% pC16c_A = parsGLT;
pC16c_B = parsTreCycle2;
% pC16c_C = parsPGM1;
% pC16c_D = parsTPS1;
% pC16c_E = parsTPS2;
% pC16c_F = parsNTH1;
%
parComb_cell = cell(ntests,1);
% 
for i = 1:ntests
%     parComb_cell{i} = pC16c_A;
    parComb_cell{i} = pC16c_B;
%     parComb_cell{i} = pC16c_C;
%     parComb_cell{i} = pC16c_D;
%     parComb_cell{i} = pC16c_E;
%     parComb_cell{i} = pC16c_F;
end
% 

% weights: 
% + w(60) -> simulated hxk.
% + w(61) -> reference parameter set.
w_simHXK = 60;
w_refY3M1 = 61;
w_simGLT = 62;
w_simTRE = 63:66;
w_simPGM1 = 63;
w_simTPS1 = 64;
w_simTPS2 = 65;
w_simNTH1 = 66;
w_simG1P = 67;
w_simT6P = 68;
w_simUDPGlc = 69;
w_simTRE2 = 63:69;
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
    warray_cell{i}(w_refY3M1) = ones;
    warray_cell{i}(w_simTRE2) = ones * warray_list(i);
end
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
% %% PARAMETER ESTIMATION: NO MS
% % 
% % for o = selectedRuns_idx
% cluster = parcluster('local');
% pool = parpool(cluster,4);
% selPars = pC16c_B;
% setup.selPars = selPars;
% parfor o = 1:ntests
% 
%     % disp
%     disp(o);
% %     fprintf('o = %f.\n', o) %disp(o)
% 
%     % find in the precreated arrays
%     selPars = parComb_cell{o};
%     warray = warray_cell{o}; 
%     saveName = names_cell{o};
% %     setup.selPars = selPars;
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
% %         save(saveName, 'xres', 'resnorm', 'residual', 'exitflag', 'warray', 't', 'selPars');
% 
%     
% end
% delete(pool)
% % disp('parallel run below')
% % % save
% % saveName = 'FF_pE4_x1.mat';
% % save(saveName,'allRes_xres','allRes_resnorm','allRes_residual','allRes_exitflag','allRes_t','allRes_warray','allRes_selPars');
% % % % % delete(pool)
% % quit


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

%% simulate results
% 
simRes1 = cell(1,length(namesHits1));
parfor i = 1:length(namesHits1)
    disp(i)
    xSel = xAll1(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF_enrichment(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
% 
       
%%
simRes1_16d_TRE = simRes1;
simRes1_16d_xAll1 = xAll1;
simRes1_16d_namesHits1 = namesHits1;
save('tempRes_16d_TRE.mat', ...
    'simRes1_16d_TRE', 'simRes1_16d_xAll1', 'simRes1_16d_namesHits1')
%%
load('tempRes_16d_TRE.mat', ...
    'simRes1_16d_TRE', 'simRes1_16d_xAll1', 'simRes1_16d_namesHits1')


%% Visualizing results: simulation
% % selResCell = simRes1_16d_TRE([1,2:19]);
% selResCell = simRes1_16d_TRE([1,9:19]);
selResCell = simRes1_16d_TRE([1,11:19]);
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;


%% Visualizing results: parameters
nPoints = length(namesHits1) - 1;
%
idxs_temp = [1 2 3 4 5, ...
            6 7 8 9 10, ...
            11 12 13 14 15, ...
            16 17, ...
            21];
pE16a_w_sim = warray_list(idxs_temp);
pE16a_w_refY3M1 = ones(1,nPoints);

% 
pE16a_error_sim_PGM1 = zeros(1,nPoints);
pE16a_error_sim_TPS1 = zeros(1,nPoints);
pE16a_error_sim_TPS2 = zeros(1,nPoints);
pE16a_error_sim_NTH1 = zeros(1,nPoints);
    % 
    pE16a_error_sim_G1P = zeros(1,nPoints);
    pE16a_error_sim_T6P = zeros(1,nPoints);
    pE16a_error_sim_UDPGlc = zeros(1,nPoints);
    %
pE16a_error_sim = zeros(1,nPoints);
pE16a_error_pars = zeros(1,nPoints);
% 
temp_time = simRes1_16d_TRE{1}.T_FF01;
temp_pgm1 = simRes1_16d_TRE{1}.V_FF01(:,17);
temp_tps1 = simRes1_16d_TRE{1}.V_FF01(:,21);
temp_tps2 = simRes1_16d_TRE{1}.V_FF01(:,19);
temp_nth1 = simRes1_16d_TRE{1}.V_FF01(:,20);
    %
    temp_g1p = simRes1_16d_TRE{1}.Y_FF01(:,21);
    temp_t6p = simRes1_16d_TRE{1}.Y_FF01(:,26);
    temp_udpglc = simRes1_16d_TRE{1}.Y_FF01(:,24);
    %
temp_timepoints = [6:5:26, 76:50:376];
temp_PGM1exp = interp1(temp_time, temp_pgm1, temp_timepoints, 'pchip');
temp_TPS1exp = interp1(temp_time, temp_tps1, temp_timepoints, 'pchip');
temp_TPS2exp = interp1(temp_time, temp_tps2, temp_timepoints, 'pchip');
temp_NTH1exp = interp1(temp_time, temp_nth1, temp_timepoints, 'pchip');
    %
    temp_G1Pexp = interp1(temp_time, temp_g1p, temp_timepoints, 'pchip');
    temp_T6Pexp = interp1(temp_time, temp_t6p, temp_timepoints, 'pchip');
    temp_UDPGlc = interp1(temp_time, temp_udpglc, temp_timepoints, 'pchip');
    %
for i = 1:nPoints
    %
    temp_PGM1sim2 = interp1(simRes1_16d_TRE{i+1}.T_FF01, simRes1_16d_TRE{i+1}.V_FF01(:,17), temp_timepoints, 'pchip');
    temp_TPS1sim2 = interp1(simRes1_16d_TRE{i+1}.T_FF01, simRes1_16d_TRE{i+1}.V_FF01(:,17), temp_timepoints, 'pchip');
    temp_TPS2sim2 = interp1(simRes1_16d_TRE{i+1}.T_FF01, simRes1_16d_TRE{i+1}.V_FF01(:,17), temp_timepoints, 'pchip');
    temp_NTH1sim2 = interp1(simRes1_16d_TRE{i+1}.T_FF01, simRes1_16d_TRE{i+1}.V_FF01(:,17), temp_timepoints, 'pchip');
        %
        temp_G1Psim2 = interp1(simRes1_16d_TRE{i+1}.T_FF01, simRes1_16d_TRE{i+1}.Y_FF01(:,21), temp_timepoints, 'pchip');
        temp_T6Psim2 = interp1(simRes1_16d_TRE{i+1}.T_FF01, simRes1_16d_TRE{i+1}.Y_FF01(:,26), temp_timepoints, 'pchip');
        temp_UDPGlcsim2 = interp1(simRes1_16d_TRE{i+1}.T_FF01, simRes1_16d_TRE{i+1}.Y_FF01(:,24), temp_timepoints, 'pchip');
        % 
    % 
    pE16a_error_sim_PGM1(i) = sum(abs( ( temp_PGM1exp - temp_PGM1sim2 ) ./ temp_PGM1exp ));
    pE16a_error_sim_TPS1(i) = sum(abs( ( temp_TPS1exp - temp_TPS1sim2 ) ./ temp_TPS1exp ));
    pE16a_error_sim_TPS2(i) = sum(abs( ( temp_TPS2exp - temp_TPS2sim2 ) ./ temp_TPS2exp ));
    pE16a_error_sim_NTH1(i) = sum(abs( ( temp_NTH1exp - temp_NTH1sim2 ) ./ temp_NTH1exp ));
        % 
        pE16a_error_sim_G1P(i) = sum(abs( ( temp_G1Pexp - temp_G1Psim2 ) ./ temp_G1Pexp ));
        pE16a_error_sim_T6P(i) = sum(abs( ( temp_T6Pexp - temp_T6Psim2 ) ./ temp_T6Pexp ));
        pE16a_error_sim_UDPGlc(i) = sum(abs( ( temp_UDPGlc - temp_UDPGlcsim2 ) ./ temp_UDPGlc ));
        %
    % 
    pE16a_error_sim(i) = pE16a_error_sim_PGM1(i) + pE16a_error_sim_TPS1(i) + pE16a_error_sim_TPS2(i) + pE16a_error_sim_NTH1(i) + ...
                            pE16a_error_sim_G1P(i) + pE16a_error_sim_T6P(i) + pE16a_error_sim_UDPGlc(i);
    % 
    pE16a_error_pars(i) = sum(abs(xAll1(i+1,[parsPGM1, parsTPS1, parsTPS2, parsNTH1]) - x_Y3M1([parsPGM1, parsTPS1, parsTPS2, parsNTH1])));
end


%% pE16c_A
% pE16c_A
figure(101)

% semilogx
subplot(2,2,1)
yyaxis left
semilogx(pE16a_w_sim, pE16a_error_sim, 'o-')
hold on
yyaxis right
semilogx(pE16a_w_sim, pE16a_error_pars, 'o-')
legend('error simTRE','error parameters')
title('w_{refY3M1} = ones -> w_{TRE} (semilogx)')
xlim([min(pE16a_w_sim) max(pE16a_w_sim)])
hold off

% loglog
subplot(2,2,2)
yyaxis left
loglog(pE16a_w_sim, pE16a_error_sim, 'o-')
hold on
yyaxis right
loglog(pE16a_w_sim, pE16a_error_pars, 'o-')
legend('error simTRE','error parameters')
title('w_{refY3M1} = ones -> w_{TRE} (loglog)')
xlim([min(pE16a_w_sim) max(pE16a_w_sim)])
hold off

% roc-like curve
subplot(2,2,3)
% area([pE16a_error_pars_18(3) pE16a_error_pars_18(2)], [10 10], ...
%     'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
hold on
plot(pE16a_error_pars, pE16a_error_sim, ...
    'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k')
set(gca,'xdir','reverse')
xlabel('error_{pars}')
ylabel('error_{simTRE}')
title('roc-like curve')
hold off

% suptitle
suptitle('regulatization for 16d, TRE')


%% Visualizing results: simulation ( close range)
% % % selResCell = simRes1_16d_TRE([1,2:19]);
% % selResCell = simRes1_16d_TRE([1,9:19]);
% selResCell = simRes1_16d_TRE([1,6:11]);
selResCell = simRes1_16d_TRE([1,8,9,10]);
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;

%%
% clf(102)
figure(102),
plot(xAll1(1,parsTreCycle2),'k.-','MarkerSize',2)
hold on
% plot(xAll1(10,parsTreCycle2),'b.-')
plot(xAll1(9,parsTreCycle2),'b.-')
plot(x_Y3M1(end,parsTreCycle2),'r.-','MarkerSize',2)
legend('Y3M2','chosen','Y3M1')
hold off


%%
% layout
maxWidth = 5;
maxRange = 1;
nSel = 3;
colRange = cool(nSel);
% xSel = xAll1([1 9 2],:);
xSel = xAll1([1 10 2],:);
% xSel = xAll1([1 19 2],:);
% xSel = xAll1([1 3 2],:);
% mBlue = [0, 0.4470, 0.7410];
% mOrange = [0.8500, 0.3250, 0.0980];

%

clf(203)
fh = figure(203);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% xlen = 1:length(parsGLT);
% tempLabel = cell(length(parsGLT),1);
% subplot(2,3,1)
% for i = xlen
%     wid = maxWidth * (abs(xSel(3,parsGLT(i)) - xSel(4,parsGLT(i)))/maxRange);
%     line([xlen(i) xlen(i)],[min([xSel(3,parsGLT(i)) xSel(4,parsGLT(i))]) max([xSel(3,parsGLT(i)) xSel(4,parsGLT(i))])],...
%         'Color','black','LineWidth',wid)
%     if(i < 3)
%         tempLabel{i} = strtok(legenda.parameters{parsGLT(i)}(10:end),', ');
%     else
%         tempLabel{i} = strtok(legenda.parameters{parsGLT(i)}(7:end),', ');
%     end
% end
% hold on
% %
% for i = 1:nSel
%     plot(xlen,abs(xSel(i,parsGLT)),'color',colRange(i,:),'LineStyle',':')
%     plot(xlen,-abs(xSel(i,parsGLT)),'color',colRange(i,:),'LineStyle',':')
% end
% %
% h = zeros(1,nSel);
% for i = 1:nSel
%     h(i) = scatter(xlen,xSel(i,parsGLT),'filled', ...
%         'MarkerEdgeColor', colRange(i,:), ...
%         'MarkerFaceColor', colRange(i,:));
% end
% xticks(xlen)
% xticklabels(tempLabel)
% ylim([-3 3])
% % legend(h,'xY3M2','x16','x15','xY3M1',...
% %     'Location','SouthOutside','Orientation','horizontal')
% title('parsGLT')
% hold off
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% xlen = 1:length(parsGLK);
% tempLabel = cell(length(parsGLK),1);
% subplot(2,3,2)
% for i = xlen
%     wid = maxWidth * (abs(xSel(3,parsGLK(i)) - xSel(4,parsGLK(i)))/maxRange);
%     line([xlen(i) xlen(i)],[min([xSel(3,parsGLK(i)) xSel(4,parsGLK(i))]) max([xSel(3,parsGLK(i)) xSel(4,parsGLK(i))])],...
%         'Color','black','LineWidth',wid)
%     tempLabel{i} = strtok(legenda.parameters{parsGLK(i)}(9:end),', ');
% end
% hold on
% %
% for i = 1:nSel
%     plot(xlen,abs(xSel(i,parsGLK)),'color',colRange(i,:),'LineStyle',':')
%     plot(xlen,-abs(xSel(i,parsGLK)),'color',colRange(i,:),'LineStyle',':')
% end
% %
% h = zeros(1,nSel);
% for i = 1:nSel
%     h(i) = scatter(xlen,xSel(i,parsGLK),'filled', ...
%         'MarkerEdgeColor', colRange(i,:), ...
%         'MarkerFaceColor', colRange(i,:));
% end
% xticks(xlen)
% xticklabels(tempLabel)
% ylim([-3 3])
% % legend(h,'xY3M2','x16','x15','xY3M1',...
% %     'Location','SouthOutside','Orientation','horizontal')
% title('parsGLK')
% hold off
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% xlen = 1:length(parsAXP);
% tempLabel = cell(length(parsAXP),1);
% subplot(2,3,3)
% for i = xlen
%     wid = maxWidth * (abs(xSel(3,parsAXP(i)) - xSel(4,parsAXP(i)))/maxRange);
%     if xSel(3,parsAXP(i)) == xSel(4,parsAXP(i))
%     else
%         line([xlen(i) xlen(i)],[min([xSel(3,parsAXP(i)) xSel(4,parsAXP(i))]) max([xSel(3,parsAXP(i)) xSel(4,parsAXP(i))])],...
%             'Color','black','LineWidth',wid)
%     end
%     tempLabel{i} = legenda.parameters{parsAXP(i)}(8:end-2);
% end
% hold on
% %
% for i = 1:nSel
%     plot(xlen,abs(xSel(i,parsAXP)),'color',colRange(i,:),'LineStyle',':')
%     plot(xlen,-abs(xSel(i,parsAXP)),'color',colRange(i,:),'LineStyle',':')
% end
% %
% h = zeros(1,nSel);
% for i = 1:nSel
%     h(i) = scatter(xlen,xSel(i,parsAXP),'filled', ...
%         'MarkerEdgeColor', colRange(i,:), ...
%         'MarkerFaceColor', colRange(i,:));
% end
% xticks(xlen)
% xticklabels(tempLabel)
% ylim([-3 3])
% % legend(h,'xY3M2','x16','x15','xY3M1',...
% %     'Location','SouthOutside','Orientation','horizontal')
% title('parsAXP')
% hold off
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
xlen = 1:length(parsTre);
tempLabel = cell(length(parsTre),1);
% subplot(2,3,[4 5])
for i = xlen
    wid = maxWidth * (abs(xSel(2,parsTre(i)) - xSel(3,parsTre(i)))/maxRange);
    line([xlen(i) xlen(i)],[min([xSel(2,parsTre(i)) xSel(3,parsTre(i))]) max([xSel(2,parsTre(i)) xSel(3,parsTre(i))])],...
        'Color','black','LineWidth',wid)
    if((i < 5)||(i == 6)||(i == 9)||(i == 12))
        tempLabel{i} = strtok(legenda.parameters{parsTre(i)}(8:end),', ');
    elseif((i == 5)||(i == 8)||(i == 10)||(i == 11))
        tempLabel{i} = legenda.parameters{parsTre(i)}(8:end-3);
    elseif((i == 7)||(i == 13)||(i == 14))
        tempLabel{i} = legenda.parameters{parsTre(i)}(8:end-2);
    end
end
hold on
%
for i = 1:nSel
    plot(xlen,abs(xSel(i,parsTre)),'color',colRange(i,:),'LineStyle',':')
    plot(xlen,-abs(xSel(i,parsTre)),'color',colRange(i,:),'LineStyle',':')
end
%
h = zeros(1,nSel);
for i = 1:nSel
    h(i) = scatter(xlen,xSel(i,parsTre),'filled', ...
        'MarkerEdgeColor', colRange(i,:), ...
        'MarkerFaceColor', colRange(i,:));
end
% edits to highlights limits
% plot(xlen,xSel(4,parsTre), 'bo:')
plot(xlen,xSel(2,parsTre), 'ko-')
xlim([0 15])
line([0 15],[-1 -1],'Color','black','LineStyle','--')
line([0 15],[1 1],'Color','black','LineStyle','--')


xticks(xlen)
xticklabels(tempLabel)
ylim([-3 3])
% legend(h,'xY3M2','x16','x15','xY3M1',...
%     'Location','SouthOutside','Orientation','horizontal')
title('parsTre')
hold off
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%

%% Visualizing results: simulation ( close range)
% % % selResCell = simRes1_16d_TRE([1,2:19]);
% % selResCell = simRes1_16d_TRE([1,9:19]);
% selResCell = simRes1_16d_TRE([1,6:11]);
selResCell = simRes1_16d_TRE([1,10]);
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;
%%

x16d_initial = xAll1(1,:);
x16d_final = xAll1(10,:);
save('x16d.mat','x16d_initial','x16d_final')

%% Starting simulation
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x16d_final,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % 
clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;
% % 
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;
