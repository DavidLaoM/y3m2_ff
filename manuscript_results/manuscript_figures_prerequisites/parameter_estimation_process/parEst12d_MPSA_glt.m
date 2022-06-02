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
plotMode = 0;

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
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
plotMode = 2; % multiple simulations
% clear selResCell
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;

plotMode = 2;
referencePlotSimulations
plotMode = 0;

safecopy_baseSim = selResCell;
clear selResCell
% 
% setup.clamp_GLCec = 0;
% 
% % % %%
% % % figure,
% % % % plot(T_FF01_1, Y_FF01_1(:,36), '.-')
% % % plot(T, Y(:,36), '.-')
% % % title('glc')


%% prepare the parameter set array

% parameter sets
parsGLT = [35 36 38 166];
% 
nMPSA = 1000;
rng(1), randVals_1 = -3 + 6 * rand(nMPSA,3);
rng(2), randVals_2 = -6 + 12 * rand(nMPSA,1);
randVals = [randVals_1, randVals_2];
% 
xMPSA = zeros(nMPSA,166);
for i = 1:nMPSA
    xMPSA(i,:) = x;
    xMPSA(i,parsGLT) = randVals(i,:);
end

% % visualize
% figure,
% subplot(1,4,1), hist(xMPSA(:,parsGLT(1))),
% subplot(1,4,2), hist(xMPSA(:,parsGLT(2))),
% subplot(1,4,3), hist(xMPSA(:,parsGLT(3))),
% subplot(1,4,4), hist(xMPSA(:,parsGLT(4))),


%% run the simulations
% save in chunks of 100 simulations
simRes1 = cell(1,nMPSA);

% 1/10
rangeTemp = 1:100;
parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_1 = simRes1(rangeTemp);
save('pE12d_simRes_1.mat','simRes1_1')

% 2/10
rangeTemp = 101:200;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_2 = simRes1(rangeTemp);
save('pE12d_simRes_2.mat','simRes1_2')

% 3/10
rangeTemp = 201:300;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_3 = simRes1(rangeTemp);
save('pE12d_simRes_3.mat','simRes1_3')

% 4/10
rangeTemp = 301:400;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_4 = simRes1(rangeTemp);
save('pE12d_simRes_4.mat','simRes1_4')

% 5/10
rangeTemp = 401:500;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_5 = simRes1(rangeTemp);
save('pE12d_simRes_5.mat','simRes1_5')

% 6/10
rangeTemp = 501:600;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_6 = simRes1(rangeTemp);
save('pE12d_simRes_6.mat','simRes1_6')

% 7/10
rangeTemp = 601:700;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_7 = simRes1(rangeTemp);
save('pE12d_simRes_7.mat','simRes1_7')

% 8/10
rangeTemp = 701:800;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_8 = simRes1(rangeTemp);
save('pE12d_simRes_8.mat','simRes1_8')

% 9/10
rangeTemp = 801:900;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_9 = simRes1(rangeTemp);
save('pE12d_simRes_9.mat','simRes1_9')

% 10/10
rangeTemp = 901:1000;
% parpool(4)
parfor i = rangeTemp%1:nMPSA
    disp(i);
    xSel = xMPSA(i,:);
    [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
    simRes1{i}.T_FF01 = T_FF01;
    simRes1{i}.Y_FF01 = Y_FF01;
    simRes1{i}.V_FF01 = V_FF01;
end
simRes1_10 = simRes1(rangeTemp);
save('pE12d_simRes_10.mat','simRes1_10')


%% laod and put all data together
% load them
load('pE12d_simRes_1.mat','simRes1_1')
load('pE12d_simRes_2.mat','simRes1_2')
load('pE12d_simRes_3.mat','simRes1_3')
load('pE12d_simRes_4.mat','simRes1_4')
load('pE12d_simRes_5.mat','simRes1_5')

load('pE12d_simRes_6.mat','simRes1_6')
load('pE12d_simRes_7.mat','simRes1_7')
% load('pE12d_simRes_8.mat','simRes1_8')
% load('pE12d_simRes_9.mat','simRes1_9')
% load('pE12d_simRes_10.mat','simRes1_10')

% put all arrays rogether
simRes1 = [simRes1_1, simRes1_2, simRes1_3, simRes1_4, simRes1_5,...
    simRes1_6, simRes1_7];
clear simRes1_1 simRes1_2 simRes1_3 simRes1_4 simRes1_5 simRes1_6 simRes1_7


% %% test simulation
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
% selResCell = simRes1_1;
% 
% referencePlotSimulations
% plotMode = 0;


%% screening for values
% vGLT at time = 20 s
% vGLT at time(end)
lenS = length(simRes1);
vGLT_20 = zeros(1,lenS);
vGLT_end = zeros(1,lenS);
for i = 1:lenS
    % Select time
    tSel = 20;
    [~,closestIndex] = min(abs(simRes1{i}.T_FF01-tSel));
    % select point
    vGLT_20(i) = real(simRes1{i}.V_FF01(closestIndex,1));
    vGLT_end(i) = real(simRes1{i}.V_FF01(end,1));
end
%%
clf(10)
figure(10)
% 
subplot(3,2,[1 2 3 4])
% range in area plot
% x = 0:15;
x = 0:0.001:0.999;
y1 = x * 0.8139/3.787;
y2 = x * 3.132/3.315;
s2 = patch([x fliplr(x)], [y1 fliplr(y2)], [.9 .9 .9]);
s2.EdgeColor = [1 1 1];

% scatter data
hold on
s = scatter(vGLT_20,vGLT_end);
s.LineWidth = 0.6;
s.SizeData = 15;
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'k';
% xlabel('v_{GLT,20}')
% ylabel('v_{GLT,end}')
% ylabel('v_{GLT,20}')
% xlabel('v_{GLT,end}')
ylabel('vHXT at 20 s (mM s^{-1})')
xlabel('vHXT at 0 or 400 s (mM s^{-1})')
xlim([0 1])%50])
ylim([0 1])

% specific data point
hold on
plot(dataset.FF01.fluxes{1}(5), dataset.FF01.fluxes{1}(2), ...
    'ko','MarkerFaceColor','b')
% legend('MPSA range', 'MPSA data points','experimental')
box on

% 
subplot(3,2,5)
hold on
for i = 1:700
    temp_x = simRes1{i}.T_FF01;
    temp_y = simRes1{i}.Y_FF01(:,36);
    plot(temp_x, temp_y, ...
        'k-', 'Linewidth', 0.5)
end
plot(safecopy_baseSim{1}.T_FF01, safecopy_baseSim{1}.Y_FF01(:,36), ...
    'b-', 'Linewidth', 1.5)
% 
s = scatter(dataset.FF01.metabolites.ECglucose.time,...
    dataset.FF01.metabolites.ECglucose.conc);
s.LineWidth = 0.6;
s.SizeData = 20;
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'r';
ylabel('GLCec concentration (mM)')
xlabel('time (s)')
ylim([0 0.6])
xlim([0 400])
box on

% 
subplot(3,2,6)
hold on
for i = 1:100
    temp_x = simRes1{i}.T_FF01;
    temp_y = simRes1{i}.V_FF01(:,1);
%     if(isreal(temp_y) == 0)
    if((isreal(temp_y) == 0)||(max(temp_y)) > 2)
        continue
    end
    plot(temp_x, temp_y, ...
        'k-', 'Linewidth', 0.5)
end
plot(safecopy_baseSim{1}.T_FF01, safecopy_baseSim{1}.V_FF01(:,1), ...
    'b-', 'Linewidth', 1.5)
% 
s = scatter(dataset.FF01.fluxes_times,...
    dataset.FF01.fluxes{1});
s.LineWidth = 0.6;
s.SizeData = 20;
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'r';
ylabel('vHXT (mM s^{-1})')
xlabel('time (s)')
ylim([0 3])
xlim([0 400])
box on

% % 
set(10, 'Position', [100 100 550 650])


%% Saving outcome
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters panel_11
print -dpdf -painters panel_2


%% memoryDump
% % x=0:0.1:10;
% x=0:15;
% % y1=exp(-x/2);
% y1 = x * 0.8139/3.787;
% % y2=exp(-x/3);
% y2 = x * 3.132/3.315;
% % clf(11)
% figure(11)
% hold all
% plot(x,y1)
% plot(x,y2)
% s2 = patch([x fliplr(x)], [y1 fliplr(y2)], [.9 .9 .9]);
% s2.EdgeColor = [1 1 1];
% % s2.LineStyle = 'none';
% hold off
% xlim([0 15])
% ylim([0 3.5])

%% memoryDump
% 
% %%
% for i = 1:16
%     fid = fopen('parEst12b_TREexport_estimationGLTGLK.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d; % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE12b_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end

