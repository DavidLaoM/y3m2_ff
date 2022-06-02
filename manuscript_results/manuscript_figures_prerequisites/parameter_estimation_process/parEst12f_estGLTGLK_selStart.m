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

setup.changing_Keq_glt = 2;
x(166) = 0; % % <== change x166 if needed here.


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
% 
% % % %%
% % % figure,
% % % % plot(T_FF01_1, Y_FF01_1(:,36), '.-')
% % % plot(T, Y(:,36), '.-')
% % % title('glc')


%% selecting the calues that showed up in parEst12e
idxs = [201 175 79 435 379]; % of 5

% parameter sets
parsGLT = [35 36 38 166];
parsGLK = 28:34;
% 
nMPSA = 1000;
rng(1), randVals_1 = -3 + 6 * rand(nMPSA,3);
rng(2), randVals_2 = -6 + 12 * rand(nMPSA,1);
randVals = [randVals_1, randVals_2];
rng(3), randVals_3 = -3 + 6 * rand(nMPSA,7);
% 
xMPSA = zeros(nMPSA,166);
for i = 1:nMPSA
    xMPSA(i,:) = x;
    xMPSA(i,parsGLT) = randVals(i,:);
    xMPSA(i,parsGLK) = randVals_3(i,:);
end

xAll_start = xMPSA(idxs,:);

%%

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
ntests = 5;


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
names_cell{1} = 'FF_pE12f_start1_nS%d.mat';
names_cell{2} = 'FF_pE12f_start2_nS%d.mat';
names_cell{3} = 'FF_pE12f_start3_nS%d.mat';
names_cell{4} = 'FF_pE12f_start4_nS%d.mat';
names_cell{5} = 'FF_pE12f_start5_nS%d.mat';
%
optimOptions2.names_cell = names_cell;


% par combinations cell: to select the parameters to optimize
parsGLT_keq = 166;
parsGLT = [35 36 38 166];
parsGLK = 28:34; % hxk
% 
pC1 = parsGLT_keq;
pC2 = parsGLT;
pC3 = parsGLK;
pC4 = [parsGLT, parsGLK];
%
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
parComb_cell{1} = pC4;
parComb_cell{2} = pC4;
parComb_cell{3} = pC4;
parComb_cell{4} = pC4;
parComb_cell{5} = pC4;
%
optimOptions2.parComb_cell = parComb_cell;
%     parComb_cell{i} = [parsHXK,parsPGM1,parsUGP,...
%         parsTPS1,parsTPS2,parsNTH1,parsGLY2];


% weights:
w_glt_glk = [39 40];
w_udp_glc = 24;
% w_glc_e = 6;
% w_g1p = 21;
% w_udpg = 24;
% w_t6p_tre = [26 25];
% w_treRates_glk_pgi = [55 59 57 58 56 40 41];
% w_gly = [86 87 88];
% w_tre_ic = 25;
% w_tre_ec = 37;
% w_ath1 = [83 84];
%     w_ath1_ec = 83;
%     w_ath1_vac = 84; 
% w_agt1 = 85;


% % % 
% % w_idxs1 = [w_tre_ic, w_tre_ec, w_ath1, w_agt1, ...
% %     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi];
% % w_idxs2 = w_idxs1; w_idxs2(end-3) = [];
% % w_idxs3 = [w_tre_ic, ...
% %     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi]; w_idxs3(end-3) = [];
% % w_idxs4 = [...
% %     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi]; w_idxs4(end-3) = [];
% % w_idxs5 = [w_g1p, w_udpg, w_gly];
% % w_idxs6 = w_gly;
% % 
% % % % w_idxs1 = w_glc_e;
% % % % w_idxs2 = [w_idxs1, w_ath1_ec, w_agt1, w_tre_ec];
% % % % w_idxs3 = [w_idxs2, w_ath1_vac, w_tre_ic];
% % % % w_idxs4 = [w_idxs3, w_treRates_glk_pgi];
% % % % w_idxs5 = [w_idxs4, w_g1p, w_udpg, w_t6p_tre];
% % % % w_idxs6 = [w_idxs5, 39];
% % w_idxs1 = w_glc_e;
% % w_idxs2 = [w_glc_e, w_tre_ec, w_treRates_glk_pgi];
% % % w_idxs1 = [w_glc_e, w_tre_ec];
% % % w_idxs2 = [w_glc_e, w_tre_ec, 39, 40];
% % % w_idxs3 = [w_glc_e, w_tre_ec, w_treRates_glk_pgi];
% w_idxs1 = w_treRates_glk_pgi;
% w_idxs2 = [w_treRates_glk_pgi, w_glc_e];
% w_idxs3 = [w_treRates_glk_pgi, w_glc_e, w_tre_ec, w_t6p_tre, w_g1p, w_udpg, 12];
% % w_idxs1 = [w_ath1_vac, w_agt1];
% % w_idxs2 = [w_ath1_vac, w_agt1, w_ath1_ec, w_tre_ic, w_tre_ec];
% % w_idxs3 = [w_ath1_vac, w_agt1, w_ath1_ec, w_tre_ic, w_tre_ec, 58];
w_idxs1 = w_glt_glk;
w_idxs2 = [w_glt_glk, w_udp_glc];
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
%

% creating blank parameters array
nS_num = 50;
rng(1), randVals = -0.25 + 0.50 * rand(nS_num,166);


%% Selected runs
% 1*16
selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS

% 1
if selCase == 1
    range_selectedRuns = 1;
    range_nS_num = 1:17;
elseif selCase == 2
    range_selectedRuns = 1;
    range_nS_num = 18:34;
elseif selCase == 3
    range_selectedRuns = 1;
    range_nS_num = 35:50;
 
% 2
elseif selCase == 4
    range_selectedRuns = 2;
    range_nS_num = 1:17;
elseif selCase == 5
    range_selectedRuns = 2;
    range_nS_num = 18:34;
elseif selCase == 6
    range_selectedRuns = 2;
    range_nS_num = 35:50;   
 
% 3
elseif selCase == 7
    range_selectedRuns = 3;
    range_nS_num = 1:17;
elseif selCase == 8
    range_selectedRuns = 3;
    range_nS_num = 18:34;
elseif selCase == 9
    range_selectedRuns = 3;
    range_nS_num = 35:50;
    
% 4
elseif selCase == 10
    range_selectedRuns = 4;
    range_nS_num = 1:17;
elseif selCase == 11
    range_selectedRuns = 4;
    range_nS_num = 18:34;
elseif selCase == 12
    range_selectedRuns = 4;
    range_nS_num = 35:50;
    
    
% 5
elseif selCase == 13
    range_selectedRuns = 5;
    range_nS_num = 1:17;
elseif selCase == 14
    range_selectedRuns = 5;
    range_nS_num = 18:34;
elseif selCase == 15
    range_selectedRuns = 5;
    range_nS_num = 35:50;
    
end

%
disp(range_selectedRuns)
disp(range_nS_num)


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
for o = range_selectedRuns %selectedRuns_idx
    for o2 = range_nS_num %1:nS_num
%         % 
%         sprintf('o = %d, o2 = %d.', o, o2)
        % 

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
%         names_cell{5} = 'FF_pE12f_start5_nS%d.mat';
        saveName = sprintf(tempName,o2);
%         saveName = sprintf(tempName,o2);
    %     saveName = [tempName, '.mat'];
%         saveName = tempName;

        % core options:
        plength = length(selPars);
        x_temp = x(selPars);
            % added randomness
            x_temp = x_temp + randVals(o2,pC4);
        
        % boundaries
        lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
        ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;
            % specific for glt
            lb(4) = -6;
            ub(4) = 0;
            
%         x_temp = x(selPars);
%     % %     lb = -3*ones(1,plength); lb(end-8:end) = -5 * ones;
%     % %     ub = 3*ones(1,plength); ub(end-8:end) = 5 * ones;
%     %     lb = -3*ones(1,plength); lb(end-12:end) = -5 * ones;
%     %     ub = 3*ones(1,plength); ub(end-12:end) = 5 * ones;
%         lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
%         ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;
%         if(o >= 13)
%             lb(end-13:end) = -10;
%             ub(end-13:end) = 10;
%         end

        % FF parEst extra
    %     NumberCycles = 20;

    %     % %% run check
    %     [error]=costfunSystemY3M1_FF_pRun(x_temp,canelas_SS,setup,xarray,data,dataset,NumberCycles,IC,selPars,warray);
    %     errorAnalysis_Y3M1;
    %     % %%

        % parameter estimation
%         sprintf('run num.%d tested',o)
        sprintf('StartCase = %d, NumRepeat = %d.', o, o2)
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
end
% disp('parallel run below')
% % save
% saveName = 'FF_pE4_x1.mat';
% save(saveName,'allRes_xres','allRes_resnorm','allRes_residual','allRes_exitflag','allRes_t','allRes_warray','allRes_selPars');
% % % % delete(pool)
% quit


% %% recall the results
% % using loop this time
% xAll1 = x;
% namesHits1 = cell(1,1); namesHits1{1} = 'initial';
% for i = 1:ntests
%     for u = 1:2
%         for j = 1:nS_num
%             %
%             loadName = sprintf(names_cell{i}, j, u);
%             %
%             if exist(loadName,'file') == 2 % if exist
%                 load(loadName);
%                 selPars = parComb_cell{i};
%                 x3 = x;
%                 x3(selPars) = xres;
%                 xAll1 = [xAll1; x3];
%                 namesHits1 = [namesHits1; loadName];
%             end
%         end
%         % naming (specific):
%         if i == 1
%             xAll1_1 = xAll1;
%             namesHits1_1 = namesHits1;
%         elseif i == 2
%             xAll1_2 = xAll1;
%             namesHits1_2 = namesHits1;
%         elseif i == 3
%             xAll1_3 = xAll1;
%             namesHits1_3 = namesHits1;
%         elseif i == 4
%             xAll1_4 = xAll1;
%             namesHits1_4 = namesHits1;
%         elseif i == 5
%             xAll1_5 = xAll1;
%             namesHits1_5 = namesHits1;
%         elseif i == 6
%             xAll1_6 = xAll1;
%             namesHits1_6 = namesHits1;
%         elseif i == 7
%             xAll1_7 = xAll1;
%             namesHits1_7 = namesHits1;
%         elseif i == 8
%             xAll1_8 = xAll1;
%             namesHits1_8 = namesHits1;
%         end
%         % restarting
%         clear xAll1 namesHits1
%         xAll1 = x;
%         namesHits1 = cell(1,1); namesHits1{1} = 'initial';
%     end
% end
% %
% xAll = [xAll1_1; xAll1_2; xAll1_3; xAll1_4; xAll1_5; xAll1_6; xAll1_7; xAll1_8]; 
% namesHits_pC1 = [namesHits1_1; namesHits1_2; namesHits1_3; namesHits1_4; namesHits1_5; namesHits1_6; namesHits1_7; namesHits1_8]; 
% 
% 
% %% simulations
% nParts = length(namesHits_pC1);
% simRes1 = cell(1,nParts);
% parpool(4)
% parfor i = 1:nParts
%     disp(i);
%     xSel = xAll(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% save('pE12c_simRes.mat','simRes1')
% 
% 
% %% visualizing the results
% % % 
% % parsGLT_keq = 166;
% % parsGLT = [35 36 38 166];
% % parsGLK = 28:34; % hxk
% % % 
% % pC1 = parsGLT_keq;
% % pC2 = parsGLT;
% % pC3 = parsGLK;
% % pC4 = [parsGLT, parsGLK];
% % load('pE12b_simRes.mat','simRes1');
% 
% % %%
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% % [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(xAll1_7(2,:),canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % % clear selResCell
% % selResCell{1}.T_FF01 = T_FF01_1;
% % selResCell{1}.Y_FF01 = Y_FF01_1;
% % selResCell{1}.V_FF01 = V_FF01_1;
% 
% % tempRange = [1:51]; % 1
% % tempRange = [52:102];
% % tempRange = [103:121];
% % tempRange = [122:139];
% % tempRange = [140:150]; % 5
% % tempRange = [151:162];
% % tempRange = [163:167];
% tempRange = [168:174];
% 
% selResCell = simRes1(tempRange);
% 
% referencePlotSimulations
% plotMode = 0;




%% memoryDump

% % %%
% for i = 1:15
%     fid = fopen('parEst12f_estGLTGLK_selStart.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d; % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE12f_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end

