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

% % % % load('x_Comb32_flx');

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

% % % % %% updated parameter sets
% % % % load('x_Comb32_flx'); x = xarray;
% % % % load('pset_20210325.mat')
% % % % load('parSet_99b.mat')
% % % % %
% % % % parsGlyco = [28:38];
% % % % parsTre = [83:86,119:128];
% % % % parsAXP = [109:111,129:131];
% % % % parsAdd = 144:158;
% % % % %
% % % % x1 = x;
% % % % x2 = x99b;
% % % %     x2(parsGlyco) = x(parsGlyco);
% % % %     x2(parsTre) = x(parsTre);
% % % % %     x2(parsAXP) = x(parsAXP);
% % % %     x2(parsAdd) = x(parsAdd);
% % % % x3 = x2;
% % % %     x3([109 129]) = x([109 129]);
% % % % clear x, x = x3;


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
% get experimental data,

% IC(40) = 100; % initial conncetration
setup.glycSynthDeg_separate = 1;

% %% sims
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% % plotflag = 2; % variables by iteration
% choosedataset
% setup.experiment = 1;
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x1,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% [T_FF01_2,Y_FF01_2,V_FF01_2] = simulate_FF(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % clear selResCell
% selResCell{1}.T_FF01 = T_FF01_1;
% selResCell{1}.Y_FF01 = Y_FF01_1;
% selResCell{1}.V_FF01 = V_FF01_1;
% selResCell{2}.T_FF01 = T_FF01_2;
% selResCell{2}.Y_FF01 = Y_FF01_2;
% selResCell{2}.V_FF01 = V_FF01_2;
% 
% referencePlotSimulations
% plotMode = 0;

%%
% save('pE7_startSims.mat','selResCell')
% load('pE7_startSims.mat','selResCell')


%% recall last case
x = x1;
x(162) = 10;
setup.ATHinhibitionT6P = 1;
load('pset_pE7_xres.mat','x_pE7_start','x_pE7_end');
clear x x1 x2
x = x_pE7_end;

% %% sims
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


% %% PSA to see if we can still keep the same params for ath1_ec and ath1_vac
% 
% nS = 16;
% xArray = zeros(nS,length(x));
% for i = 1:nS
%     xArray(i,:) = x;
%     if i <= 8 % kcat increasing
%         xArray(i,150) = xArray(i,150) + (i-1);
%     else % km decreasing
%         xArray(i,149) = xArray(i,149) - (i-9);
%     end
% end
% 
% % simulations in parfor loop
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% plotMode = 0;
% 
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % pool = parpool(cluster,4);
% parfor o = 1:nS
%     % simulate
%     x_sel = xArray(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE8_psa_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
%     
% end
% delete(pool)
% quit
% 
% 
% %% plot PSA
% simRes1 = cell(1,nS);
% for i = 1:nS
%     loadName = sprintf('pE8_psa_n%d.mat',i);
%     load(loadName);
%     simRes1(i) = selResCell;
% end
% % simRes1(1) = [];
% simRes1_a = simRes1(1:8);
% simRes1_b = simRes1(9:end);
% 
% 
% %% plot data
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% 
% choosedataset
% 
% % plotMode = 1; % single simulation
% plotMode = 2; % multiple simulations
% % plotMode = 10; % all
% 
% % simRes1_a
% selResCell = simRes1_a;
% referencePlotSimulations
% % simRes1_b
% selResCell = simRes1_b;
% referencePlotSimulations


%% separating ath1_ec and ath1_vac
x(163:164) = x(149:150); % km_tre, k_cat
x(165) = x(162); % ki_t6p

    % 149_Ktre_ATHvac->163_Ktre_ATHec
    % 150_kcat_ATHvac->164_kcat_ATHec
setup.ATH_separate_EC_VAC = 1;


% %% sims
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
% % disp('something')


%% PSA over the new reactions added
x_interval = -5:1:5;
nrand = 11*3;
rng(1), x_rand = 5 - 10 * rand(nrand,3);

xPSA = zeros(length(x_interval)*3+nrand,length(x));
[len_PSA,~] = size(xPSA);
for i = 1:len_PSA
    xPSA(i,:) = x;
    if((i >= 1)&&(i <= 11))
        i2 = i - 0;
        xPSA(i,163) = x_interval(i2);
    elseif((i >= 12)&&(i <= 22))
        i2 = i - 11;
        xPSA(i,164) = x_interval(i2);
    elseif((i >= 23)&&(i <= 33))
        i2 = i - 22;
        xPSA(i,165) = x_interval(i2);
    else
        i2 = i - 33;
        xPSA(i,163:165) = xPSA(i2,163:165);
    end
end
%% simulations in parfor loop
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
plotMode = 0;

% cluster = parcluster('local');
% % pool = parpool(cluster,16);
% pool = parpool(cluster,4);
% parfor o = 1:len_PSA
%     disp(o)
%     % simulate
%     x_sel = xPSA(o,:);
%     [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x_sel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     selResCell = cell(1,1);
%     selResCell{1}.T_FF01 = T_FF01_1;
%     selResCell{1}.Y_FF01 = Y_FF01_1;
%     selResCell{1}.V_FF01 = V_FF01_1;
%     % save
%     saveName = sprintf('pE8_psa_2_n%d.mat',o);
%     parsave_Y3M2_cluster_2(saveName, selResCell);
%     
% end
% delete(pool)
% quit


% %% plot PSA
% simRes1 = cell(1,len_PSA);
% for i = 1:len_PSA
%     loadName = sprintf('pE8_psa_2_n%d.mat',i);
%     load(loadName);
%     simRes1(i) = selResCell;
% end
% % simRes1(1) = [];
% simRes1_a = simRes1(1:11);
% simRes1_b = simRes1(12:22);
% simRes1_c = simRes1(23:33);
% simRes1_d = simRes1(34:end);
% 
% 
% %% plot data
% [legenda] = legendaFull; %legenda for the names needed
% metNames = legenda.metabolites;
% reactNames = legenda.fluxes;
% 
% choosedataset
% 
% % plotMode = 1; % single simulation
% plotMode = 2; % multiple simulations
% % plotMode = 10; % all
% 
% % simRes1_a
% selResCell = simRes1_a;
% referencePlotSimulations
% % simRes1_b
% selResCell = simRes1_b;
% referencePlotSimulations
% % simRes1_c
% selResCell = simRes1_c;
% referencePlotSimulations
% % simRes1_d
% selResCell = simRes1_d;
% referencePlotSimulations


%% Parameter estimation setup
setup.fakeATHcostFuns = 1;

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
ntests = 30;

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
    names_cell{i} = sprintf('FF_pE9b_wA%d',i);
end
optimOptions2.names_cell = names_cell;

% par combinations cell: to select the parameters to optimize
parsHXK = 28:34; % hxk
parsPGM1 = 83:86; % pgm1
parsUGP = 144:148; % ugp
parsTPS1 = 124:128; % tps1
parsTPS2 = 119:121; % tps2
parsNTH1 = 122:123; % nth1
% parsGLY = 159;
parsGLY2 = 160:161;
%     parsNTH1 = 122:123;       % nth1      [122:123]
%     parsATH1 = 149:150;       % ath1_v.e  [149:150]
    parsATH1 = [149:150,162:165];       % ath1_v.e  [149:150]
        parsATH1vac = [149:150,162]; % ath1_v.vac  [149:150]
        parsATH1ec = 163:165; % ath1_v.ec  [149:150]
    parsAGT1 = [151:154,158];   % agt1      [151:154,158]
    parsVACT = 155:156;       % vacT      [155:156]

parCombs1 = parsATH1vac;
parCombs2 = [parsATH1vac, parsATH1ec];
parCombs3 = [parsATH1vac, parsAGT1];
parCombs4 = [parsATH1vac, parsVACT];
parCombs5 = [parsATH1vac, parsATH1ec, parsAGT1, parsVACT];
    
% 
parComb_cell = cell(ntests,1);
% for i = 1:ntests
%     if i <= 2
%         parComb_cell{i} = [parsHXK,parsPGM1,parsUGP,...
%             parsTPS1,parsTPS2,parsNTH1,parsGLY2,...
%             parsATH1,parsAGT1,parsVACT];
%     elseif i <= 17
%         parComb_cell{i} = [parsNTH1,parsATH1,parsAGT1,parsVACT];
%     else
%         parComb_cell{i} = parsAGT1;
%     end
% end
for i = 1:ntests
%     if i <= 2
%         parComb_cell{i} = [parsHXK,parsPGM1,parsUGP,...
%             parsTPS1,parsTPS2,parsNTH1,parsGLY2,...
%             parsATH1,parsAGT1,parsVACT];
%     elseif i <= 16
    if i <= 16
        parComb_cell{i} = [parsNTH1,parsATH1,parsAGT1,parsVACT];
    else
        parComb_cell{i} = parsATH1;
    end
end
% parCombs1
parComb_cell{1} = parCombs1;
parComb_cell{2} = parCombs1;
parComb_cell{3} = parCombs1;
parComb_cell{4} = parCombs1;
parComb_cell{5} = parCombs1;
parComb_cell{6} = parCombs1;
% parCombs2
parComb_cell{7} = parCombs2;
parComb_cell{8} = parCombs2;
parComb_cell{9} = parCombs2;
parComb_cell{10} = parCombs2;
parComb_cell{11} = parCombs2;
parComb_cell{12} = parCombs2;
% parCombs3
parComb_cell{13} = parCombs3;
parComb_cell{14} = parCombs3;
parComb_cell{15} = parCombs3;
parComb_cell{16} = parCombs3;
parComb_cell{17} = parCombs3;
parComb_cell{18} = parCombs3;
% parCombs4
parComb_cell{19} = parCombs4;
parComb_cell{20} = parCombs4;
parComb_cell{21} = parCombs4;
parComb_cell{22} = parCombs4;
parComb_cell{23} = parCombs4;
parComb_cell{24} = parCombs4;
% parCombs5
parComb_cell{25} = parCombs5;
parComb_cell{26} = parCombs5;
parComb_cell{27} = parCombs5;
parComb_cell{28} = parCombs5;
parComb_cell{29} = parCombs5;
parComb_cell{30} = parCombs5;
% 
optimOptions2.parComb_cell = parComb_cell;
%     parComb_cell{i} = [parsHXK,parsPGM1,parsUGP,...
%         parsTPS1,parsTPS2,parsNTH1,parsGLY2];

% weights:
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
% 
w_idxs1 = [w_tre_ic, w_tre_ec, w_ath1, w_agt1, ...
    w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi];
w_idxs2 = w_idxs1; w_idxs2(end-3) = [];
% w_idxs3 = [w_tre_ic, w_tre_ec, ...
%     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi];
% w_idxs3(end-3) = [];
% 
% w_idxs4 = [w_tre_ic, ...
%     w_g1p, w_udpg, w_t6p_tre, w_gly, w_treRates_glk_pgi];
% w_idxs4(end-3) = [];

w_array1 = w_ath1_ec;
w_array2 = [w_ath1_ec, w_ath1_vac];
w_array3 = [w_ath1_ec, w_ath1_vac, w_tre_ec];
w_array4 = [w_ath1_ec, w_ath1_vac, w_tre_ec, w_tre_ic];
w_array5 = [w_ath1_ec, w_ath1_vac, w_tre_ec, w_tre_ic, w_treRates_glk_pgi];
w_array6 = [w_ath1_ec, w_treRates_glk_pgi];


% %
% rng(1), randArr = 1 + (10-1) * rand(ntests,length(w_idxs2));
% 
warray_cell = cell(ntests,1);
% base
% for i = 1:ntests % base,
%     warray_cell{i} = blankWeight;
%     if((i == 1))
%         warray_cell{i}(w_idxs1) = ones;
%     elseif((i == 2))
%         warray_cell{i}(w_idxs2) = ones;
%     elseif((i == 3) || (i == 18))
%         warray_cell{i}(w_idxs2) = ones;
%     else
%         warray_cell{i}(w_idxs2) = randArr(i,:); 
%     end
% end
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
% % %     
% %     if((i == 1))
% %         warray_cell{i}(w_idxs1) = ones;
% %     elseif((i == 2))
% %         warray_cell{i}(w_idxs2) = ones;
% %     elseif((i == 3) || (i == 18))
% %         warray_cell{i}(w_idxs2) = ones;
% %     else
% %         warray_cell{i}(w_idxs2) = randArr(i,:); 
% %     end
end

% parCombs1
warray_cell{1}(w_array1) = ones;
warray_cell{2}(w_array2) = ones;
warray_cell{3}(w_array3) = ones;
warray_cell{4}(w_array4) = ones;
warray_cell{5}(w_array5) = ones;
warray_cell{6}(w_array6) = ones;
% parCombs2
warray_cell{7}(w_array1) = ones;
warray_cell{8}(w_array2) = ones;
warray_cell{9}(w_array3) = ones;
warray_cell{10}(w_array4) = ones;
warray_cell{11}(w_array5) = ones;
warray_cell{12}(w_array6) = ones;
% parCombs3
warray_cell{13}(w_array1) = ones;
warray_cell{14}(w_array2) = ones;
warray_cell{15}(w_array3) = ones;
warray_cell{16}(w_array4) = ones;
warray_cell{17}(w_array5) = ones;
warray_cell{18}(w_array6) = ones;
% parCombs4
warray_cell{19}(w_array1) = ones;
warray_cell{20}(w_array2) = ones;
warray_cell{21}(w_array3) = ones;
warray_cell{22}(w_array4) = ones;
warray_cell{23}(w_array5) = ones;
warray_cell{24}(w_array6) = ones;
% parCombs5
warray_cell{25}(w_array1) = ones;
warray_cell{26}(w_array2) = ones;
warray_cell{27}(w_array3) = ones;
warray_cell{28}(w_array4) = ones;
warray_cell{29}(w_array5) = ones;
warray_cell{30}(w_array6) = ones;


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
    saveName = [tempName, '.mat'];
    
    % core options:
    x_temp = x(selPars);
    plength = length(selPars);
% %     lb = -3*ones(1,plength); lb(end-8:end) = -5 * ones;
% %     ub = 3*ones(1,plength); ub(end-8:end) = 5 * ones;
%     lb = -3*ones(1,plength); lb(end-12:end) = -5 * ones;
%     ub = 3*ones(1,plength); ub(end-12:end) = 5 * ones;
    lb = -5*ones(1,plength); %lb(end-12:end) = -5 * ones;
    ub = 5*ones(1,plength); %ub(end-12:end) = 5 * ones;
    
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
    loadName = [names_cell{i}, '.mat'];
    if exist(loadName,'file') == 2 % if exist
        load(loadName);
        selPars = parComb_cell{i};
        x3 = x;
        x3(selPars) = xres;
        xAll1 = [xAll1; x3];
        namesHits1 = [namesHits1; loadName];
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
save('pE9b_simRes.mat','simRes1')
% save('pE8_simRes_updated.mat','simRes1')
%%
nParts = length(namesHits1);
simRes1 = cell(1,nParts);
load('pE9b_simRes.mat','simRes1')


%% plot data
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;

choosedataset

% plotMode = 1; % single simulation
plotMode = 2; % multiple simulations
% plotMode = 10; % all

% simRes_pC1
selResCell = simRes1;
referencePlotSimulations

% %%
% % simRes_pC1: first part. [parsNTH1,parsATH1,parsAGT1,parsVACT];
% selResCell = simRes1([1,2:17]);
% referencePlotSimulations
% % simRes_pC1: second part. parsATH1
% selResCell = simRes1([1,18:33]);
% referencePlotSimulations


% % % % %%
% % % % for i = 1:length(simRes1)
% % % %     TREvac120 = interp1(selResCell{i}.T_FF01, selResCell{i}.Y_FF01(:,38), 120, 'pchip');
% % % %     if TREvac120 > 25
% % % %         fprintf('i = %f, TREvac120 = %f.\n', i, TREvac120);
% % % %     end
% % % % %     UDPglc200 = interp1(selResCell{i}.T_FF01, selResCell{i}.Y_FF01(:,24), 200, 'pchip');
% % % % %     if UDPglc200 > 1.6
% % % % %         fprintf('i = %f, UDPglc200 = %f.\n', i, UDPglc200);
% % % % %     end
% % % % end
% % % % 
% % % % 
% % % % %%
% % % % % select
% % % % selIdxs = [4,5,12,13];
% % % % % selIdxs = 18;
% % % % selIdxs2 = 1:length(simRes1); selIdxs2(selIdxs) = [];
% % % % % 1
% % % % tempSave_lab = namesHits1;
% % % % namesHits1 = namesHits1(selIdxs);
% % % % selResCell = simRes1(selIdxs);
% % % % % plot
% % % % referencePlotSimulations
% % % % % recall
% % % % selResCell = simRes1;
% % % % namesHits1 = tempSave_lab;
% % % % % 2
% % % % tempSave_lab = namesHits1;
% % % % namesHits1 = namesHits1(selIdxs2);
% % % % selResCell = simRes1(selIdxs2);
% % % % % plot
% % % % referencePlotSimulations
% % % % % recall
% % % % selResCell = simRes1;
% % % % namesHits1 = tempSave_lab;
% % % % 
% % % % 
% % % % %% more focuse selection, based on the baseline of Tre_vac
% % % % for i = selIdxs
% % % %     TREvac120 = interp1(selResCell{i}.T_FF01, selResCell{i}.Y_FF01(:,38), 120, 'pchip');
% % % % %     if TREvac120 > 25
% % % %         fprintf('i = %f, TREvac120 = %f.\n', i, TREvac120);
% % % % %     end
% % % % end
% % % % 
% % % % 
% % % % %%
% % % % % select
% % % % selIdxs3 = 12;
% % % % % 1
% % % % tempSave_lab = namesHits1;
% % % % namesHits1 = namesHits1(selIdxs3);
% % % % selResCell = simRes1(selIdxs3);
% % % % % plot
% % % % referencePlotSimulations
% % % % % recall
% % % % selResCell = simRes1;
% % % % namesHits1 = tempSave_lab;
% % % % 
% % % % 
% % % % %% probably the case is num.6, w.7
% % % % x_pE6_x1_start = xAll1(1,:);
% % % % x_pE6_x1_end = xAll1(12,:);
% % % % % 
% % % % save('pset_pE6_x1res.mat','x_pE6_x1_start','x_pE6_x1_end');


% % %%
% for i = 1:30
%     fid = fopen('parEst9b_ath_ec_more_focused.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d; % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE9b_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end
