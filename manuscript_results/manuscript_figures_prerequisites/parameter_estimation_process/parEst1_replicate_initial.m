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
plotflag = 0;
% plotflag = 2; % variables by iteration
choosedataset
setup.experiment = 1;
% [T_FF01,Y_FF01,V_FF01] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotflag = 0;


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
    'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
%     'OutputFcn',{@saveIterationsMain},...
%     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% setup.parEst.costfun = 10015; 

nParCombs = 1;
nWeightCombs = 10;
nExtras = 0;
nMS = 0;
ntests = nParCombs * nWeightCombs + nExtras;

% names cell (to save)
names_cell = cell(ntests,1);
tempName = 'FF_pE_pCs%d_wA%d'; %.mat to be added upon call
for i = 1:ntests
    pC_num = fix((i-1)/nWeightCombs) + 1;
    warr_num = rem(i,nWeightCombs);
    if warr_num == 0, warr_num = nWeightCombs; end
    names_cell{i} = sprintf(tempName,pC_num,warr_num);
end
optimOptions2.names_cell = names_cell;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% par combinations cell: to select the parameters to optimize
parsHXK = 28:34;
parComb_cell = cell(ntests,1);
parComb_cell{1} = parsHXK; %28:34;
parComb_cell{2} = parsHXK; %28:34;
parComb_cell{3} = parsHXK; %28:34;
parComb_cell{4} = parsHXK; %28:34;
parComb_cell{5} = parsHXK; %28:34;
parComb_cell{6} = parsHXK; %28:34;
parComb_cell{7} = parsHXK; %28:34;
parComb_cell{8} = parsHXK; %28:34;
parComb_cell{9} = parsHXK; %28:34;
parComb_cell{10} = parsHXK; %28:34;
%
optimOptions2.parComb_cell = parComb_cell;

% weights:
w_hxk = 40;
w_pep = 12;
    
warray_cell = cell(ntests,1);
% base
for i = 1:ntests % base,
    warray_cell{i} = blankWeight;
    warray_cell{i}(w_hxk) = ones;
 end
% specific
warray_cell{1}(w_pep) = 0;
warray_cell{2}(w_pep) = 1E-2;
warray_cell{3}(w_pep) = 3E-2;
warray_cell{4}(w_pep) = 1E-1;
warray_cell{5}(w_pep) = 3E-1;

warray_cell{6}(w_pep) = 1E0;
warray_cell{7}(w_pep) = 3E0;
warray_cell{8}(w_pep) = 1E1;
warray_cell{9}(w_pep) = 3E1;
warray_cell{10}(w_pep) = 1E2;


%% Selected runs
% 1*16
selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS
selectedRuns_idx = selCase
% selectedRuns_idx = [1:6] + 6 * (selCase - 1)


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


%% DATA RECALL (not the multistart cases yet)
% xFullScale_99b_parCombs5_warray8.mat is lacking at start
xAll1 = x;
namesHits1 = cell(1,1); namesHits1{1} = 'initial';
for i = 1:ntests
    loadName = [names_cell{i}, '.mat'];
    if exist(loadName,'file') == 2 % if exist
        load(loadName);
        selPars = parComb_cell{i};
        x2 = x;
        x2(selPars) = xres;
        xAll1 = [xAll1; x2];
        namesHits1 = [namesHits1; loadName];
    end
end


%% simulate results
nParts = length(namesHits1);
simRes1 = cell(1,nParts);
for i = 1:nParts
    disp(i);
    xSel = xAll1(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end


%% plot data
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;

choosedataset

% plotMode = 1; % single simulation
plotMode = 2; % multiple simulations
% plotMode = 10; % all

selResCell = simRes1;
referencePlotSimulations


% %%
% % 16 * 6
% % separate in different parCombs
% xAll1_pC1 = xAll1([1,(1:6)+1],:);
% xAll1_pC2 = xAll1([1,(1:6)+7],:);
% xAll1_pC3 = xAll1([1,(1:5)+13],:);
% 
% % names
% namesHits1_pC1 = namesHits1([1,(1:6)+1]);
% namesHits1_pC2 = namesHits1([1,(1:6)+7]);
% namesHits1_pC3 = namesHits1([1,(1:5)+13]);


% %% Recall params and labels + simulations
% % %% individual + gp only
% 
% disp(1);
% [simRes_pC1] = simulateY3M1_separateGPSS(xAll1_pC1, canelas_SS, data, dataset, setup);




%%
% % % %%
% for i = 1:10
%     fid = fopen('parEst1_replicate_initial.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE1_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end


