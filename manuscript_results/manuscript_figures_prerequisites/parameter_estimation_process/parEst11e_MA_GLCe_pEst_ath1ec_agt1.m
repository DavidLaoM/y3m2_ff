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
% IC(40) = 100; % initial conncetration
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

    % 149_Ktre_ATHvac->163_Ktre_ATHec
    % 150_kcat_ATHvac->164_kcat_ATHec
setup.ATH_separate_EC_VAC = 1;


%% setting the option
setup.glycogenReactionsSink = 1;
setup.dataset = dataset;
%
load('pset_pE10_xres.mat','x_pE10_start','x_pE10_end'); x = x_pE10_end;

%% latest setup
setup.updated_bmf_Cx_ATH1ec = 1;
setup.updated_bmf_Cx_ATH1ec2 = 1;
x(166) = x(162);


% %% leaving the clamp on ATHec_AGT1 ON
% % 
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
% % selResCell(1) = []
% referencePlotSimulations
% plotMode = 0;
% 
% 
% %% NTH1 OFF, see if tre_ic accumulater and others behave 
% % tic
% x2 = x; x2(123) = -10;
% %
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{2}.T_FF01 = T_FF01_1;
% selResCell{2}.Y_FF01 = Y_FF01_1;
% selResCell{2}.V_FF01 = V_FF01_1;
% % selResCell(1) = []
% referencePlotSimulations
% plotMode = 0;
% % toc
% 
% 
% %% NTH1 OFF, see if tre_ic accumulater and others behave 
% % tic
% x3 = x2; x3(164) = x3(164)+5; % % % % x2(123) = -10;
% %
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x3,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{3}.T_FF01 = T_FF01_1;
% selResCell{3}.Y_FF01 = Y_FF01_1;
% selResCell{3}.V_FF01 = V_FF01_1;
% % selResCell(1) = []
% referencePlotSimulations
% plotMode = 0;


%% Parameter estimation setup

% blank and constant setup
blankWeight = zeros(1,88); % 85+3 for glycerol
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
names_cell{1} = 'FF_pE11e_pC1_wA1.mat';
names_cell{2} = 'FF_pE11e_pC2_wA1.mat';
names_cell{3} = 'FF_pE11e_pC3_wA1.mat';
names_cell{4} = 'FF_pE11e_pC4_wA1.mat';
names_cell{5} = 'FF_pE11e_pC5_wA1.mat';
names_cell{6} = 'FF_pE11e_pC6_wA1.mat';
names_cell{7} = 'FF_pE11e_pC7_wA1.mat';
names_cell{8} = 'FF_pE11e_pC8_wA1.mat';
%
names_cell{9} = 'FF_pE11e_pC1_wA2.mat';
names_cell{10} = 'FF_pE11e_pC2_wA2.mat';
names_cell{11} = 'FF_pE11e_pC3_wA2.mat';
names_cell{12} = 'FF_pE11e_pC4_wA2.mat';
names_cell{13} = 'FF_pE11e_pC5_wA2.mat';
names_cell{14} = 'FF_pE11e_pC6_wA2.mat';
names_cell{15} = 'FF_pE11e_pC7_wA2.mat';
names_cell{16} = 'FF_pE11e_pC8_wA2.mat';
%
optimOptions2.names_cell = names_cell;

% par combinations cell: to select the parameters to optimize
parsHXK = 28:34; % hxk
parsPGM1 = 83:86; % pgm1
parsUGP = 144:148; % ugp
parsTPS1 = 124:128; % tps1
parsTPS2 = 119:121; % tps2
parsNTH1 = 122:123; % nth1
    parsTreCycle = [parsPGM1, parsUGP, parsTPS1, parsTPS2, parsNTH1];
% parsGLY = 159;
parsGLY2 = 160:161;
%     parsNTH1 = 122:123;       % nth1      [122:123]
%     parsATH1 = 149:150;       % ath1_v.e  [149:150]
    parsATH1 = [149:150,162:165];       % ath1_v.e  [149:150]
        parsATH1vac = [149 150 162];
        parsATH1ec = [163 164 165];
    parsAGT1 = [151:154,158,166];   % agt1      [151:154,158]
    parsVACT = 155:156;       % vacT      [155:156]
parsGLT = [35 36 38];
    
pC1 = parsATH1ec;
pC2 = [pC1, parsAGT1];
pC3 = [pC2, parsATH1vac, parsVACT];
pC4 = [pC3, parsGLT];
pC5 = [pC3, parsHXK];
pC6 = [pC3, parsTreCycle];
pC7 = [pC3, parsGLT, parsHXK];
pC8 = [pC3, parsGLT, parsHXK, parsTreCycle];
% parCombs1 = [parsATH1ec, parsAGT1];
% % parCombs2 = [parCombs1, parsATH1vac, parsVACT];
% % parCombs3 = [parCombs2, parsPGM1, parsUGP, parsTPS1, parsTPS2, parsNTH1];
% % parCombs4 = [parCombs3, parsHXK];
% % parCombs5 = [parCombs4, 35 36 38];
% parCombs2 = [parCombs1, 35 36 38];
% parCombs3 = [parCombs1, parsHXK];
% parCombs4 = [parCombs1, 35 36 38, parsHXK];
% parCombs5 = [parCombs1, 35 36 38, parsHXK, parsPGM1, parsUGP, parsTPS1, parsTPS2, parsNTH1];
    
parComb_cell = cell(ntests,1);
% for i = 1:ntests
% %     if i <= 2
% %         parComb_cell{i} = [parsHXK,parsPGM1,parsUGP,...
% %             parsTPS1,parsTPS2,parsNTH1,parsGLY2,...
% %             parsATH1,parsAGT1,parsVACT];
% %     elseif i <= 16
%     if i <= 16
%         parComb_cell{i} = [parsNTH1,parsATH1,parsAGT1,parsVACT];
%     else
%         parComb_cell{i} = parsATH1;
%     end
% end
% 
% parComb_all = [parsHXK,parsPGM1,parsUGP,...
%     parsTPS1,parsTPS2,parsNTH1,parsGLY2,...
%     parsATH1,parsAGT1,parsVACT];
% parComb_GlyTre = [parsHXK,parsPGM1,parsUGP,...
%     parsTPS1,parsTPS2,parsNTH1,parsGLY2];
% parComb_GlyOnly = [parsPGM1,parsUGP];
% 
parComb_cell{1} = pC1;
parComb_cell{2} = pC2;
parComb_cell{3} = pC3;
parComb_cell{4} = pC4;
parComb_cell{5} = pC5;
parComb_cell{6} = pC6;
parComb_cell{7} = pC7;
parComb_cell{8} = pC8;
%
parComb_cell{9} = pC1;
parComb_cell{10} = pC2;
parComb_cell{11} = pC3;
parComb_cell{12} = pC4;
parComb_cell{13} = pC5;
parComb_cell{14} = pC6;
parComb_cell{15} = pC7;
parComb_cell{16} = pC8;
%
optimOptions2.parComb_cell = parComb_cell;
%     parComb_cell{i} = [parsHXK,parsPGM1,parsUGP,...
%         parsTPS1,parsTPS2,parsNTH1,parsGLY2];

% weights:
w_glc_e = 6;
w_g1p = 21;
w_udpg = 24;
w_t6p_tre = [26 25];
w_treRates_glk_pgi = [55 59 57 58 56 40 41];
w_gly = [86 87 88];
w_tre_ic = 25;
w_tre_ec = 37;
w_ath1 = [83 84];
    w_ath1_ec = 83;
    w_ath1_vac = 84; 
w_agt1 = 85;
% % 
% w_idxs1 = [w_tre_ic, w_tre_ec, w_ath1, w_agt1, ...
%     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi];
% w_idxs2 = w_idxs1; w_idxs2(end-3) = [];
% w_idxs3 = [w_tre_ic, ...
%     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi]; w_idxs3(end-3) = [];
% w_idxs4 = [...
%     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi]; w_idxs4(end-3) = [];
% w_idxs5 = [w_g1p, w_udpg, w_gly];
% w_idxs6 = w_gly;
% 
% % w_idxs1 = w_glc_e;
% % w_idxs2 = [w_idxs1, w_ath1_ec, w_agt1, w_tre_ec];
% % w_idxs3 = [w_idxs2, w_ath1_vac, w_tre_ic];
% % w_idxs4 = [w_idxs3, w_treRates_glk_pgi];
% % w_idxs5 = [w_idxs4, w_g1p, w_udpg, w_t6p_tre];
% % w_idxs6 = [w_idxs5, 39];
w_idxs1 = w_glc_e;
w_idxs2 = [w_glc_e, w_tre_ec, w_treRates_glk_pgi];
% w_idxs1 = [w_glc_e, w_tre_ec];
% w_idxs2 = [w_glc_e, w_tre_ec, 39, 40];
% w_idxs3 = [w_glc_e, w_tre_ec, w_treRates_glk_pgi];


% w_idxs1 = [w_ath1_vac, w_agt1];
% w_idxs2 = [w_ath1_vac, w_agt1, w_ath1_ec, w_tre_ic, w_tre_ec];
% w_idxs3 = [w_ath1_vac, w_agt1, w_ath1_ec, w_tre_ic, w_tre_ec, 58];
% 
% %
% rng(1), randArr = 1 + (10-1) * rand(ntests,length(w_idxs2));
% 
warray_cell = cell(ntests,1);
% % base
% % for i = 1:ntests % base,
% %     warray_cell{i} = blankWeight;
% %     if((i == 1))
% %         warray_cell{i}(w_idxs1) = ones;
% %     elseif((i == 2))
% %         warray_cell{i}(w_idxs2) = ones;
% %     elseif((i == 3) || (i == 18))
% %         warray_cell{i}(w_idxs2) = ones;
% %     else
% %         warray_cell{i}(w_idxs2) = randArr(i,:); 
% %     end
% % end
for i = 1:ntests % base,
    warray_cell{i} = blankWeight;
%     if((i == 1) || (i == 17))
%         warray_cell{i}(83) = ones;
%     elseif((i == 2) || (i == 18))
%         warray_cell{i}(w_ath1) = ones;
%     elseif((i == 2) || (i == 18))
%         warray_cell{i}(w_idxs2) = ones;
%     else
%         warray_cell{i}(w_idxs2) = randArr(i,:); 
%     end
end
% 
warray_cell{1}(w_idxs1) = ones;
warray_cell{2}(w_idxs1) = ones;
warray_cell{3}(w_idxs1) = ones;
warray_cell{4}(w_idxs1) = ones;
warray_cell{5}(w_idxs1) = ones;
warray_cell{6}(w_idxs1) = ones;
warray_cell{7}(w_idxs1) = ones;
warray_cell{8}(w_idxs1) = ones;
%
warray_cell{9}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
warray_cell{10}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
warray_cell{11}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
warray_cell{12}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
warray_cell{13}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
warray_cell{14}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
warray_cell{15}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
warray_cell{16}(w_idxs2) = [1  1/5  1/10 1/10 1/10 1/10 1/10 1/10 1/10];
%


%% Selected runs
% 1*16
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
% % allRes_xres = cell(1,ntests);
% % allRes_resnorm = allRes_xres;
% % allRes_residual = allRes_xres;
% % allRes_exitflag = allRes_xres;
% % allRes_warray = allRes_xres;
% % allRes_t = allRes_xres;
% % allRes_selPars = allRes_xres;
% global blankValue
% global history
% % % % cluster = parcluster('local');
% % % % pool = parpool(cluster,16);
% % % % parfor o = selectedRuns_idx
for o = selectedRuns_idx
    disp(o)
%     blankValue = [];
%     history.x = [];%.x; x];
%     history.iteration = [];%.iteration; optimValues.iteration];
%     history.funccount = [];%.funccount; optimValues.funccount];
%     history.stepsize = [];%.stepsize; optimValues.stepsize];
%     history.gradient = [];%.gradient; optimValues.gradient'];
%     history.firstorderopt = [];%.firstorderopt; optimValues.firstorderopt];
%     history.cgiterations = [];%.cgiteration; optimValues.cgiteration];
%     history.positivedefinite = [];%.positivedefinite; optimValues.positivedefinite];
%     history.ratio = [];%.ratio; optimValues.ratio];
%     history.degenerate = [];%.degenerate; optimValues.degenerate];
%     history.trustregionradius = [];%.trustregionradius; optimValues.trustregionradius];
%     history.residual = [];%.residual; optimValues.residual'];
%     history.resnorm = [];%.resnorm; optimValues.resnorm];
       
    % find in the precreated arrays
    selPars = parComb_cell{o};
%     setup.caseStudy.parameters = selPars;
    warray = warray_cell{o}; 
%     setup.w = warray;
    tempName = names_cell{o};
%     saveName = [tempName, '.mat'];
    saveName = tempName;
    
    % core options:
    x_temp = x(selPars);
    plength = length(selPars);
% %     lb = -3*ones(1,plength); lb(end-8:end) = -5 * ones;
% %     ub = 3*ones(1,plength); ub(end-8:end) = 5 * ones;
%     lb = -3*ones(1,plength); lb(end-12:end) = -5 * ones;
%     ub = 3*ones(1,plength); ub(end-12:end) = 5 * ones;
    lb = -10*ones(1,plength); %lb(end-12:end) = -5 * ones;
    ub = 10*ones(1,plength); %ub(end-12:end) = 5 * ones;
    if((o <= 3)||( (o >= 9) && (o <= 11) ))
        lb(15:end) = -3;
        ub(15:end) = 3;
    end
    
    % FF parEst extra
    NumberCycles = 20;
    
%     % %% run check
%     [error]=costfunSystemY3M1_FF_pRun(x_temp,canelas_SS,setup,xarray,data,dataset,NumberCycles,IC,selPars,warray);
%     errorAnalysis_Y3M1;
%     % %%

    % parameter estimation
    sprintf('run num.%d tested',o)
    tic
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF_pRun,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,NumberCycles,IC,selPars,warray);
    t = toc; 
    disp(xres);

%     % recall info back in cells
%     allRes_xres{o}.xres = xres;
%     allRes_resnorm{o} = resnorm;
%     allRes_residual{o} = residual;
%     allRes_exitflag{o} = exitflag;
%     allRes_warray{o} = warray;
%     allRes_t{o} = t;
%     allRes_selPars{o} = selPars;
    
    % saving
% % % %     parsave_Y3M2_cluster(saveName, xres, resnorm, residual, exitflag, warray, t, selPars);
    save(saveName, 'xres', 'resnorm', 'residual', 'exitflag', 'warray', 't', 'selPars');
end
% disp('parallel run below')
% % save
% saveName = 'FF_pE4_x1.mat';
% save(saveName,'allRes_xres','allRes_resnorm','allRes_residual','allRes_exitflag','allRes_t','allRes_warray','allRes_selPars');
% % % % delete(pool)
% quit


%% DATA RECALL (not the multistart cases yet)
% xFullScale_99b_parCombs5_warray8.mat is lacking at start
xAll1 = x;
namesHits1 = cell(1,1); namesHits1{1} = 'initial';
for i = 1:ntests
    loadName = names_cell{i};
    if exist(loadName,'file') == 2 % if exist
% % % %             loadName = [folder, '\workContinuation\tempResults\', loadName];
        load(loadName);
        selPars = parComb_cell{i};
        x3 = x;
        x3(selPars) = xres;
        xAll1 = [xAll1; x3];
        namesHits1 = [namesHits1; loadName];
% % % %             namesHits1 = [namesHits1; loadName(end-21:end)];
        
    end
end


%% simulate results
    % 1-6
    % 7-14
nParts = length(namesHits1);
simRes1 = cell(1,nParts);
parpool(4)
parfor i = 1:nParts
    disp(i);
    xSel = xAll1(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
save('pE11e_simRes.mat','simRes1')
% % quit
% % save('pE8_simRes_updated.mat','simRes1')
%%
nParts = length(namesHits1);
simRes1 = cell(1,nParts);
load('pE11e_simRes.mat','simRes1')


%% plot data
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;

choosedataset

% plotMode = 1; % single simulation
plotMode = 2; % multiple simulations
% plotMode = 10; % all

% tempSave for the entire thing
tempSave = simRes1;

% wA1
selResCell = tempSave([1,2:9]);
referencePlotSimulations

% wA2
selResCell = tempSave([1,10:17]);
referencePlotSimulations

% % all sims
% % simRes1 = tempSave;
% selResCell = simRes1;
% referencePlotSimulations

% % %%
% for i = 1:16
%     fid = fopen('parEst11e_MA_GLCe_pEst_ath1ec_agt1.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d; % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE11e_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end
