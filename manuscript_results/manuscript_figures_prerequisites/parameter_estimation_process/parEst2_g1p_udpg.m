% % PAREST1_REPLICATE_INITIAL.m
% In this code the tre cycle and glt, hxk are optimized as it happen before
% to obtain the parameter set x32.

% structure is:
    % early sampling
    % parameter estimation with paraemter and cost function combinations
    

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
clear x, x = x3;


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
plotflag = 0

% %% sims
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01,Y_FF01,V_FF01] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{1}.T_FF01 = T_FF01;
% selResCell{1}.Y_FF01 = Y_FF01;
% selResCell{1}.V_FF01 = V_FF01;
% referencePlotSimulations
% plotMode = 0;


% % % % %% PSA on hxk to see why it does not have an effect on the estimation this time
% % % % nPSA = 5;
% % % % changeVals = [0; -1; -0.5; 0.5; 1];
% % % % xAll28 = [x; x; x; x; x]; 
% % % % xAll29 = xAll28; xAll30 = xAll28; xAll31 = xAll28; 
% % % % xAll32 = xAll28; xAll33 = xAll28; xAll34 = xAll28;
% % % % 
% % % % xAll28(:,28) = xAll28(:,28) + changeVals;
% % % % xAll29(:,29) = xAll29(:,29) + changeVals;
% % % % xAll30(:,30) = xAll30(:,30) + changeVals;
% % % % xAll31(:,31) = xAll31(:,31) + changeVals;
% % % % xAll32(:,32) = xAll32(:,32) + changeVals;
% % % % xAll33(:,33) = xAll33(:,33) + changeVals;
% % % % xAll34(:,34) = xAll34(:,34) + changeVals;
% % % % 
% % % % % nParts = length(namesHits1);
% % % % simRes28 = cell(1,nPSA);
% % % % simRes29 = simRes28; simRes30 = simRes28; simRes31 = simRes28;
% % % % simRes32 = simRes28; simRes33 = simRes28; simRes34 = simRes28;
% % % % %
% % % % parfor i = 1:nPSA
% % % %     disp(i);
% % % %     % 28
% % % %     xSel = xAll28(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes28{i}.T_FF01 = T_FF01;
% % % %     simRes28{i}.Y_FF01 = Y_FF01;
% % % %     simRes28{i}.V_FF01 = V_FF01;
% % % %     % 29
% % % %     xSel = xAll29(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes29{i}.T_FF01 = T_FF01;
% % % %     simRes29{i}.Y_FF01 = Y_FF01;
% % % %     simRes29{i}.V_FF01 = V_FF01;
% % % %     % 30
% % % %     xSel = xAll30(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes30{i}.T_FF01 = T_FF01;
% % % %     simRes30{i}.Y_FF01 = Y_FF01;
% % % %     simRes30{i}.V_FF01 = V_FF01;
% % % %     % 31
% % % %     xSel = xAll31(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes31{i}.T_FF01 = T_FF01;
% % % %     simRes31{i}.Y_FF01 = Y_FF01;
% % % %     simRes31{i}.V_FF01 = V_FF01;
% % % %     % 32
% % % %     xSel = xAll32(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes32{i}.T_FF01 = T_FF01;
% % % %     simRes32{i}.Y_FF01 = Y_FF01;
% % % %     simRes32{i}.V_FF01 = V_FF01;
% % % %     % 33
% % % %     xSel = xAll33(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes33{i}.T_FF01 = T_FF01;
% % % %     simRes33{i}.Y_FF01 = Y_FF01;
% % % %     simRes33{i}.V_FF01 = V_FF01;
% % % %     % 34
% % % %     xSel = xAll34(i,:);
% % % %     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % %     simRes34{i}.T_FF01 = T_FF01;
% % % %     simRes34{i}.Y_FF01 = Y_FF01;
% % % %     simRes34{i}.V_FF01 = V_FF01;
% % % % end
% % % % %%
% % % % save('tempResSave_20210611.mat', ...
% % % %     'simRes28','simRes29','simRes30',...
% % % %     'simRes31','simRes32','simRes33',...
% % % %     'simRes34')
% % % % 
% % % % %%
% % % % load('tempResSave_20210611.mat')
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
% % % % % 28
% % % % selResCell = simRes28;
% % % % referencePlotSimulations
% % % % % 29
% % % % selResCell = simRes29;
% % % % referencePlotSimulations
% % % % % 30
% % % % selResCell = simRes30;
% % % % referencePlotSimulations
% % % % % 31
% % % % selResCell = simRes31;
% % % % referencePlotSimulations
% % % % % 32
% % % % selResCell = simRes32;
% % % % referencePlotSimulations
% % % % % 33
% % % % selResCell = simRes33;
% % % % referencePlotSimulations
% % % % % 34
% % % % selResCell = simRes34;
% % % % referencePlotSimulations


%% Parameter estimation setup

% blank and constant setup
blankWeight = zeros(1,85);
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

% optimOptions2
enzNum = 6;
nchoosek_cell = cell(enzNum,1);
for i = 1:enzNum
    nchoosek_cell{i} = nchoosek(1:enzNum,i);
end
% concatenate
if exist('z','var')
    clear z
end
for i = 1:(enzNum-1)
    if exist('z','var')
    else
        z = nchoosek_cell{i}';
    end
    B = nchoosek_cell{i+1}';
    sA = size(z);
    sB = size(B);
    C = max(sA(1),sB(1));
    z = [[z;zeros(abs([C 0]-sA))],[B;zeros(abs([C,0]-sB))]];
end
z = z';
% z = z([1:3,5:6,8,11,15],:);

% nParCombs = 1;
[nParCombs,~] = size(z);
nWeightCombs = 8; % 3 specific + 5 growing part
nExtras = 0;
nMS = 0;
ntests = nParCombs * nWeightCombs + nExtras;

% names cell (to save)
names_cell = cell(ntests,1);
tempName = 'FF_pE2_pCs%d_wA%d'; %.mat to be added upon call
for i = 1:ntests
    pC_num = fix((i-1)/nWeightCombs) + 1;
    warr_num = rem(i,nWeightCombs);
    if warr_num == 0, warr_num = nWeightCombs; end
    names_cell{i} = sprintf(tempName,pC_num,warr_num);
end
optimOptions2.names_cell = names_cell;

% par combinations cell: to select the parameters to optimize
parComb_cell = cell(ntests,1);
% parameter options
parOptions = cell(enzNum,1);
parOptions{1} = 28:34; % hxk
parOptions{2} = 83:86; % pgm1
parOptions{3} = 144:148; % ugp
parOptions{4} = 124:128; % tps1
parOptions{5} = 119:121; % tps2
parOptions{6} = 122:123; % nth1
for i = 1:ntests
    pselected = [];
    pC_num = fix((i-1)/nWeightCombs) + 1;
    for j = 1:enzNum
%         findVal = find(z(pC_num,:) == j);
        findVal = ismember(j,z(pC_num,:));
        if findVal == 1
            pselected = [pselected, parOptions{j}];
        end
    end
    parComb_cell{i} = pselected;
end
optimOptions2.parComb_cell = parComb_cell;

% weights:
w_g1p = 21;
w_udpg = 24;
w_treRates_glk_pgi = [55 59 57 58 56 40 41];
w_t6p_tre = [26 25];

warray_cell = cell(ntests,1);
% base
for i = 1:ntests % base,
    warray_cell{i} = blankWeight;
    warray_cell{i}([w_g1p w_udpg]) = ones;
end
% empty first two cases
for i = 1:nParCombs
    i1 = i*8-7;
    warray_cell{i1}(w_udpg) = zeros;
    i2 = i*8-6;
    warray_cell{i2}(w_g1p) = zeros;
end
% add the latest part
for i = 1:nParCombs
    %
    i4 = i*8-4;
    warray_cell{i4}(w_treRates_glk_pgi) = 1E-3 * ones/2/length(w_treRates_glk_pgi);
    warray_cell{i4}(w_t6p_tre) = 1E-3 * ones/2/length(w_t6p_tre);
    %
    i5 = i*8-3;
    warray_cell{i5}(w_treRates_glk_pgi) = 1E-2 * ones/2/length(w_treRates_glk_pgi);
    warray_cell{i5}(w_t6p_tre) = 1E-2 * ones/2/length(w_t6p_tre);
    %
    i6 = i*8-2;
    warray_cell{i6}(w_treRates_glk_pgi) = 1E-1 * ones/2/length(w_treRates_glk_pgi);
    warray_cell{i6}(w_t6p_tre) = 1E-1 * ones/2/length(w_t6p_tre);
    %
    i7 = i*8-1;
    warray_cell{i7}(w_treRates_glk_pgi) = 3E-1 * ones/2/length(w_treRates_glk_pgi);
    warray_cell{i7}(w_t6p_tre) = 3E-1 * ones/2/length(w_t6p_tre);
    %
    i8 = i*8-0;
    warray_cell{i8}(w_treRates_glk_pgi) = 1E0 * ones/2/length(w_treRates_glk_pgi);
    warray_cell{i8}(w_t6p_tre) = 1E0 * ones/2/length(w_t6p_tre);
end


%% Selected runs
% 1*16
selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS
% selectedRuns_idx = selCase
if selCase <= 24
    selectedRuns_idx = [1:16] + 16 * (selCase - 1);
else
    selectedRuns_idx = [1:15] + 15 * (selCase - 25) + 384;
end
disp(selectedRuns_idx)

%% PARAMETER ESTIMATION: NO MS
% global blankValue
% global history
for o = selectedRuns_idx
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
       
%     find in the precreated arrays
    selPars = parComb_cell{o};
    setup.caseStudy.parameters = selPars;
    tempName = names_cell{o};
    saveName = [tempName, '.mat'];
    warray = warray_cell{o}; 
    setup.w = warray;
        
%     core options:
    x_temp = x(selPars);
    plength = length(selPars);
    lb = -3*ones(1,plength);
    ub = 3*ones(1,plength);
    
    % FF parEst extra
    NumberCycles = 20;
    
% %     % %% run check
%     [error]=costfunSystemY3M1_FF(x_temp,canelas_SS,setup,xarray,data,dataset,NumberCycles,IC);
% %     errorAnalysis_Y3M1;
%     % %%

%     parameter estimation
    sprintf('run num.%d tested',o)
    tic
%     [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_temp_separateParsGPSS,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,warray);
%                                                function [error]=costfunSystemY3M1_temp_separateParsGPSS(x_temp,canelas_SS,setup,x,data,dataset,w)
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_FF,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,NumberCycles,IC);
%                                                function [error]=costfunSystemY3M1_FF(x_temp,canelas_SS,setup,xarray,data,dataset,n,IC0)

    t = toc; 
    disp(xres);
%     save
%     save(saveName,'xres','resnorm','residual','exitflag','t','warray','lam','history','blankValue');
    save(saveName,'xres','resnorm','residual','exitflag','t','warray');
    
end
% quit

%% DATA RECALL (not the multistart cases yet)
% xFullScale_99b_parCombs5_warray8.mat is lacking at start
xAll1 = x;
namesHits1 = cell(1,1); namesHits1{1} = 'initial';
timeSpent = zeros(1,1);
for i = 1:ntests
    loadName = [names_cell{i}, '.mat'];
    if exist(loadName,'file') == 2 % if exist
        load(loadName);
        selPars = parComb_cell{i};
        x2 = x;
        x2(selPars) = xres;
        xAll1 = [xAll1; x2];
        namesHits1 = [namesHits1; loadName];
        timeSpent = [timeSpent; t];
    end
end


%% simulate results
    % 1-6
    % 7-14
nParts = length(namesHits1);
% simRes1 = cell(1,nParts);
% parpool(3)
% parfor i=1:3, c(:,i) = eig(rand(1000)); end
% parfor i = 1:nParts
% simRes1 = cell(1,50); % simulating first 6 parameter combinations
simRes1 = cell(1,100); % simulating first 6 parameter combinations
% xAll2 = xAll1([1, 50:(50+48)],:);
% xAll2 = xAll1([1, 99:198],:);
xAll2 = xAll1([1, 199:297],:);
parfor i = 1:100
%     if i == 1, i2 = 1; else, i2 = 1 -48; end
    disp(i);
%     xSel = xAll1(i,:);
    xSel = xAll2(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
%     simRes1{i2}.T_FF01 = T_FF01;
%     simRes1{i2}.Y_FF01 = Y_FF01;
%     simRes1{i2}.V_FF01 = V_FF01;
end
%% save 1-49
% % % simRes_1_49 = simRes1;
% % % save('simRes_2_1_49.mat', 'simRes_1_49')
% % simRes_50_98 = simRes1;
% % save('simRes_2_50_98.mat', 'simRes_50_98')
% simRes_101_198 = simRes1;
% save('simRes_2_101_198.mat', 'simRes_101_198')
simRes_199_297 = simRes1;
save('simRes_2_199_297.mat', 'simRes_199_297')


%% load all of them back
load('simRes_2_1_49.mat'), temp1 = simRes_1_49;
load('simRes_2_50_98.mat'), temp2 = simRes_50_98;
load('simRes_2_101_198.mat'), temp3 = simRes_101_198;
load('simRes_2_199_297.mat'), temp4 = simRes_199_297;
simRes_actual = [temp1, temp2(2:end), temp3(2:end), temp4(2:end)];
% namesHits1_actual = namesHits1(1:98);
clear temp1 temp2 temp3 temp4
%% split them

% simRes
simRes_pC1 = simRes_actual([1,2:9]);
simRes_pC2 = simRes_actual([1,10:17]);
simRes_pC3 = simRes_actual([1,18:25]);
simRes_pC4 = simRes_actual([1,26:33]);
simRes_pC5 = simRes_actual([1,34:41]);
simRes_pC6 = simRes_actual([1,42:49]);
simRes_pC7 = simRes_actual([1,50:57]);
simRes_pC8 = simRes_actual([1,58:65]);
simRes_pC9 = simRes_actual([1,66:73]);
simRes_pC10 = simRes_actual([1,74:81]);
simRes_pC11 = simRes_actual([1,82]);
% % % simRes_pC12 = simRes_actual([1,2:9]);
simRes_pC13 = simRes_actual([1,83:90]);
simRes_pC14 = simRes_actual([1,91:98]);


simRes_pC15 = simRes_actual([1,99:106]);
simRes_pC16 = simRes_actual([1,107:114]);
simRes_pC17 = simRes_actual([1,115:122]);
simRes_pC18 = simRes_actual([1,123:130]);
simRes_pC19 = simRes_actual([1,131:138]);
simRes_pC20 = simRes_actual([1,139:145]);

simRes_pC21 = simRes_actual([1,146:153]);
simRes_pC22 = simRes_actual([1,154:161]);
simRes_pC23 = simRes_actual([1,162]);
% % simRes_pC24 = simRes_actual([1,10:17]);
simRes_pC25 = simRes_actual([1,163:170]);

simRes_pC26 = simRes_actual([1,171:173]);
simRes_pC27 = simRes_actual([1,174:181]);
simRes_pC28 = simRes_actual([1,182]);
simRes_pC29 = simRes_actual([1,183:188]);
% % simRes_pC30 = simRes_actual([1,10:17]);

simRes_pC31 = simRes_actual([1,189]);
% % simRes_pC32 = simRes_actual([1,10:17]);
simRes_pC33 = simRes_actual([1,190:197]);
simRes_pC34 = simRes_actual([1,198:199]);
simRes_pC35 = simRes_actual([1,200:207]);

simRes_pC36 = simRes_actual([1,208:210]);
simRes_pC37 = simRes_actual([1,211:218]);
simRes_pC38 = simRes_actual([1,219:224]);
simRes_pC39 = simRes_actual([1,225:232]);
simRes_pC40 = simRes_actual([1,233:238]);

simRes_pC41 = simRes_actual([1,239]);
% % simRes_pC42 = simRes_actual([1,200:207]);
simRes_pC43 = simRes_actual([1,240:247]);
% % simRes_pC44 = simRes_actual([1,200:207]);
simRes_pC45 = simRes_actual([1,248:255]);
% 
% % simRes_pC46 = simRes_actual([1,256:263]);
simRes_pC47 = simRes_actual([1,256:263]);
% % simRes_pC48 = simRes_actual([1,200:207]);
simRes_pC49 = simRes_actual([1,264]);
simRes_pC50 = simRes_actual([1,265]);

simRes_pC51 = simRes_actual([1,266]);
simRes_pC52 = simRes_actual([1,267:268]);
simRes_pC53 = simRes_actual([1,269:270]);
simRes_pC54 = simRes_actual([1,271:273]);
simRes_pC55 = simRes_actual([1,274:281]);

simRes_pC56 = simRes_actual([1,282:286]);
simRes_pC57 = simRes_actual([1,287:289]);
% % simRes_pC58 = simRes_actual([1,200:207]);
% % simRes_pC59 = simRes_actual([1,200:207]);
simRes_pC60 = simRes_actual([1,290:291]);

% % simRes_pC61 = simRes_actual([1,200:207]);
simRes_pC62 = simRes_actual([1,292:297]);
% % simRes_pC63 = simRes_actual([1,200:207]);
    % 

% loadNames
namesHits1_pC1 = namesHits1_actual([1,2:9]);
namesHits1_pC2 = namesHits1_actual([1,10:17]);
namesHits1_pC3 = namesHits1_actual([1,18:25]);
namesHits1_pC4 = namesHits1_actual([1,26:33]);
namesHits1_pC5 = namesHits1_actual([1,34:41]);
namesHits1_pC6 = namesHits1_actual([1,42:49]);
namesHits1_pC7 = namesHits1_actual([1,50:57]);
namesHits1_pC8 = namesHits1_actual([1,58:65]);
namesHits1_pC9 = namesHits1_actual([1,66:73]);
namesHits1_pC10 = namesHits1_actual([1,74:81]);
namesHits1_pC11 = namesHits1_actual([1,82]);
% % namesHits1_pC12 = namesHits1_actual([1,2:9]);
namesHits1_pC13 = namesHits1_actual([1,83:90]);
namesHits1_pC14 = namesHits1_actual([1,91:98]);



%% plot data
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;

choosedataset

% plotMode = 1; % single simulation
plotMode = 2; % multiple simulations
% plotMode = 10; % all

% selResCell = simRes1;

% simRes_pC1
selResCell = simRes_pC1;
referencePlotSimulations

% simRes_pC2
selResCell = simRes_pC2;
referencePlotSimulations

% simRes_pC3
selResCell = simRes_pC3;
referencePlotSimulations

% simRes_pC4
selResCell = simRes_pC4;
referencePlotSimulations

% simRes_pC5
selResCell = simRes_pC5;
referencePlotSimulations

% simRes_pC6
selResCell = simRes_pC6;
referencePlotSimulations


% simRes_pC7
selResCell = simRes_pC7;
referencePlotSimulations

% simRes_pC8
selResCell = simRes_pC8;
referencePlotSimulations

% simRes_pC9
selResCell = simRes_pC9;
referencePlotSimulations

% simRes_pC10
selResCell = simRes_pC10;
referencePlotSimulations

% simRes_pC11
selResCell = simRes_pC11;
referencePlotSimulations

% % simRes_pC12
% selResCell = simRes_pC12;
% referencePlotSimulations

% simRes_pC13
selResCell = simRes_pC13;
referencePlotSimulations

% simRes_pC14
selResCell = simRes_pC14;
referencePlotSimulations

%% 15-20
selResCell = simRes_pC15;
referencePlotSimulations

selResCell = simRes_pC16;
referencePlotSimulations

selResCell = simRes_pC17;
referencePlotSimulations

selResCell = simRes_pC18;
referencePlotSimulations

selResCell = simRes_pC19;
referencePlotSimulations

selResCell = simRes_pC20;
referencePlotSimulations

%% 21-25
selResCell = simRes_pC21;
referencePlotSimulations

selResCell = simRes_pC22;
referencePlotSimulations

selResCell = simRes_pC23;
referencePlotSimulations

% % selResCell = simRes_pC24;
% % referencePlotSimulations

selResCell = simRes_pC25;
referencePlotSimulations

%% 26-30
selResCell = simRes_pC26;
referencePlotSimulations

selResCell = simRes_pC27;
referencePlotSimulations

selResCell = simRes_pC28;
referencePlotSimulations

selResCell = simRes_pC29;
referencePlotSimulations

% % selResCell = simRes_pC30;
% % referencePlotSimulations

%% 31-35
selResCell = simRes_pC31;
referencePlotSimulations

% % selResCell = simRes_pC32;
% % referencePlotSimulations

selResCell = simRes_pC33;
referencePlotSimulations

selResCell = simRes_pC34;
referencePlotSimulations

selResCell = simRes_pC35;
referencePlotSimulations

%% 36-40
selResCell = simRes_pC36;
referencePlotSimulations

selResCell = simRes_pC37;
referencePlotSimulations

selResCell = simRes_pC38;
referencePlotSimulations

selResCell = simRes_pC39;
referencePlotSimulations

selResCell = simRes_pC40;
referencePlotSimulations

%% 41-45
selResCell = simRes_pC41;
referencePlotSimulations

% % selResCell = simRes_pC42;
% % referencePlotSimulations

selResCell = simRes_pC43;
referencePlotSimulations

% % selResCell = simRes_pC44;
% % referencePlotSimulations

selResCell = simRes_pC45;
referencePlotSimulations

%% 46-50
% % selResCell = simRes_pC46;
% % referencePlotSimulations

selResCell = simRes_pC47;
referencePlotSimulations

% % selResCell = simRes_pC48;
% % referencePlotSimulations

selResCell = simRes_pC49;
referencePlotSimulations

selResCell = simRes_pC50;
referencePlotSimulations

%% 51-55
selResCell = simRes_pC51;
referencePlotSimulations

selResCell = simRes_pC52;
referencePlotSimulations

selResCell = simRes_pC53;
referencePlotSimulations

selResCell = simRes_pC54;
referencePlotSimulations

selResCell = simRes_pC55;
referencePlotSimulations

%% 56-60
selResCell = simRes_pC56;
referencePlotSimulations

selResCell = simRes_pC57;
referencePlotSimulations

% % selResCell = simRes_pC58;
% % referencePlotSimulations
% % 
% % selResCell = simRes_pC59;
% % referencePlotSimulations

selResCell = simRes_pC60;
referencePlotSimulations

%% 61-63
% selResCell = simRes_pC61;
% referencePlotSimulations

selResCell = simRes_pC62;
referencePlotSimulations

% selResCell = simRes_pC63;
% referencePlotSimulations

%% recheck th interestin case is the last one
selResCell = simRes_pC62([1,end]);
referencePlotSimulations
%%


%%


nParts = length(namesHits1);
% simRes1 = cell(1,nParts);
% parpool(3)
% parfor i=1:3, c(:,i) = eig(rand(1000)); end
% parfor i = 1:nParts
% simRes1 = cell(1,50); % simulating first 6 parameter combinations
simRes1 = cell(1,100); % simulating first 6 parameter combinations
% xAll2 = xAll1([1, 50:(50+48)],:);
% xAll2 = xAll1([1, 99:198],:);
xAll2 = xAll1([1, 199:297],:);
parfor i = 1:100
%     if i == 1, i2 = 1; else, i2 = 1 -48; end
    disp(i);
%     xSel = xAll1(i,:);
    xSel = xAll2(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
%     simRes1{i2}.T_FF01 = T_FF01;
%     simRes1{i2}.Y_FF01 = Y_FF01;
%     simRes1{i2}.V_FF01 = V_FF01;
end


%%
x_pE2_start = xAll1(1,:);
x_pE2_end = xAll1(end,:);

save('pset_pE2res.mat','x_pE2_start','x_pE2_end');

% %%
% size(xAll1)



% %%
% % % %%
% for i = 1:32
%     fid = fopen('parEst2_g1p_udpg.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE2_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end


