% % 
% 

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
% % % % % load('pset_pE2res.mat');
% % % % % x1 = x_pE2_start;
% % % % % x2 = x_pE2_end;
%% glycerol_SynthDeg active
setup.glycSynthDeg = 1; % this adds a reaction to the mass balance (ODEs), reaction calculation in rate equations and the starting parameter value
% % % % % x1(159) = -4; % -10;
% % % % % x2(159) = -4; %-10;
%% extension glycogen metabolism
% load parameters
% % % % % load('pset_pE6_x1res.mat','x_pE6_x1_end'); x1 = x_pE6_x1_end;
% % % % % load('pset_pE6_x2res.mat','x_pE6_x2_end'); x2 = x_pE6_x2_end;
%% adding a more detailed glycogen metabolism
setup.glycSynthDeg_separate = 1;
%% recall last case
% % % % % x = x1;
% % % % % x(162) = 10;
setup.ATHinhibitionT6P = 1;
% % % % % load('pset_pE7_xres.mat','x_pE7_start','x_pE7_end');
% % % % % clear x x1 x2
% % % % % x = x_pE7_end;
%% separating ath1_ec and ath1_vac
% % % % % x(163:164) = x(149:150); % km_tre, k_cat
% % % % % x(165) = x(162); % ki_t6p
setup.ATH_separate_EC_VAC = 1;
%% setting the option
setup.glycogenReactionsSink = 1;
setup.dataset = dataset;
%
% % % % % load('pset_pE10_xres.mat','x_pE10_start','x_pE10_end'); x = x_pE10_end;
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
% load('x16c_E_TPS2.mat', 'x16c_E_start', 'x16c_E_final')
% x = x16c_E_final;
load('pset_Y3M2.mat','x16d_initial','x16d_final')
x = x16d_final;

%% Directly implementing the Csmin idea in out model
% general clamping and setup options
NumberCycles = 5;
% setup.clamp_GLCec = 0;
setup.clamp_GLCec = 1;
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

% 
simRes_initial = selResCell;
clear selResCell

setup.csmin = 0;


%% PART 1. MPSA

% parameter sets
parsGLT = [35 36 38];
% 
nMPSA = 1000;
rng(1), randVals = -3 + 6 * rand(nMPSA,3);
% 
xMPSA = zeros(nMPSA,length(x));
for i = 1:nMPSA
    xMPSA(i,:) = x;
    xMPSA(i,parsGLT) = randVals(i,:);
end

% simRes1 = cell(1,nMPSA);
% 
% % 1/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 1:100;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_1 = simRes1(rangeTemp);
% save('fig2_simRes_1.mat','simRes1_1')
% % 
% delete(pool)
% 
% % 2/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 101:200;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_2 = simRes1(rangeTemp);
% save('fig2_simRes_2.mat','simRes1_2')
% % 
% delete(pool)
% 
% % 3/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 201:300;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_3 = simRes1(rangeTemp);
% save('fig2_simRes_3.mat','simRes1_3')
% % 
% delete(pool)
% 
% % 4/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 301:400;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_4 = simRes1(rangeTemp);
% save('fig2_simRes_4.mat','simRes1_4')
% % 
% delete(pool)
% 
% % 5/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 401:500;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_5 = simRes1(rangeTemp);
% save('fig2_simRes_5.mat','simRes1_5')
% % 
% delete(pool)
% 
% % 6/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 501:600;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_6 = simRes1(rangeTemp);
% save('fig2_simRes_6.mat','simRes1_6')
% % 
% delete(pool)
% 
% % 7/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 601:700;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_7 = simRes1(rangeTemp);
% save('fig2_simRes_7.mat','simRes1_7')
% % 
% delete(pool)
% 
% % 8/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 701:800;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_8 = simRes1(rangeTemp);
% save('fig2_simRes_8.mat','simRes1_8')
% % 
% delete(pool)
% 
% % 9/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 801:900;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_9 = simRes1(rangeTemp);
% save('fig2_simRes_9.mat','simRes1_9')
% % 
% delete(pool)
% 
% % 10/10
% cluster = parcluster('local');
% pool = parpool(cluster,16);
% % 
% rangeTemp = 901:1000;
% % parpool(4)
% parfor i = rangeTemp%1:nMPSA
%     disp(i);
%     xSel = xMPSA(i,:);
%     [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%     simRes1{i}.T_FF01 = T_FF01;
%     simRes1{i}.Y_FF01 = Y_FF01;
%     simRes1{i}.V_FF01 = V_FF01;
% end
% simRes1_10 = simRes1(rangeTemp);
% save('fig2_simRes_10.mat','simRes1_10')
% % 
% delete(pool)


%% load and put all data together
setup.csmin = 0.094;
% load them
load('fig2_simRes_1.mat','simRes1_1')
load('fig2_simRes_2.mat','simRes1_2')
load('fig2_simRes_3.mat','simRes1_3')
load('fig2_simRes_4.mat','simRes1_4')
load('fig2_simRes_5.mat','simRes1_5')
% 
load('fig2_simRes_6.mat','simRes1_6')
load('fig2_simRes_7.mat','simRes1_7')
load('fig2_simRes_8.mat','simRes1_8')
load('fig2_simRes_9.mat','simRes1_9')
load('fig2_simRes_10.mat','simRes1_10')

% put all arrays rogether
simRes1 = [simRes1_1, simRes1_2, simRes1_3, simRes1_4, simRes1_5,...
    simRes1_6, simRes1_7, simRes1_8, simRes1_9, simRes1_10];
clear simRes1_1 simRes1_2 simRes1_3 simRes1_4 simRes1_5 
clear simRes1_6 simRes1_7 simRes1_8 simRes1_9 simRes1_10 

%
% fig2_simRes_4


%% screening for values
% 
exp_time = dataset.FF01.metabolites.ECglucose.time(1:15);
exp_GLCec = dataset.FF01.metabolites.ECglucose.conc(1:15);
% 
exp_v_time = dataset.FF01.fluxes_times;
exp_GLT = dataset.FF01.fluxes{1};
% vGLT at time = 20 s
% vGLT at time(end)
lenS = length(simRes1);
% lenS = 5;
vGLT_20 = zeros(1,lenS);
vGLT_end = zeros(1,lenS);
GLCec_end = zeros(1,lenS);
errorGLTend = zeros(1,lenS);
errorGLCend = zeros(1,lenS);
errorGLC = zeros(1,lenS);
errorGLT = zeros(1,lenS);
for i = 1:lenS
    % Select time
    tSel = 20;
    [~,closestIndex] = min(abs(simRes1{i}.T_FF01-tSel));
    % select point
    vGLT_20(i) = real(simRes1{i}.V_FF01(closestIndex,1));
    vGLT_end(i) = real(simRes1{i}.V_FF01(end,1));
    GLCec_end(i) = real(simRes1{i}.Y_FF01(end,36));
    errorGLTend(i) = abs(vGLT_end(i) - dataset.FF01.fluxes{1}(5)); %real(simRes1{i}.Y_FF01(end,36));
    errorGLCend(i) = abs(GLCec_end(i) - dataset.FF01.metabolites.ECglucose.conc(1)); %real(simRes1{i}.V_FF01(end,1));
    % 
    sim_GLCec = interp1(simRes1{i}.T_FF01, simRes1{i}.Y_FF01(:,36), exp_time,'pchip');
    errorGLC(i) = sum(abs(sim_GLCec - exp_GLCec));
    % 
    sim_GLT = interp1(simRes1{i}.T_FF01, simRes1{i}.V_FF01(:,1), exp_v_time,'pchip');
    errorGLT(i) = sum(abs(sim_GLT - exp_GLT));
end
%%
% clf(10)
figure(10)
% 
% % % % subplot(3,2,[1 2 3 4])

% % range in area plot
% % x = 0:15;
% x = 0:0.001:0.995;
% y1 = x * 0.8139/3.787;
% y2 = x * 3.132/3.315;
% s2 = patch([x fliplr(x)], [y1 fliplr(y2)], [.9 .9 .9]);
% s2.EdgeColor = [1 1 1];

% scatter data
hold on

s = scatter(vGLT_20,vGLT_end);
% s = scatter(vGLT_end, vGLT_20);
    % parsGLT = [35 36 38];
%     s = scatter(vGLT_20,vGLT_end,20, xMPSA(:,38),'Filled');

% s = scatter(vGLT_end, GLCec_end);
% s = scatter(errorGLCend, errorGLTend);
% s = scatter(errorGLC, vGLT_end);
% s = scatter(errorGLC, errorGLT);
s.LineWidth = 0.6;
s.SizeData = 15;
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'k';
% % % xlabel('error_{GLCec}')
% % % ylabel('error_{GLT}')
% % xlabel('error_{GLCec}')
% % ylabel('v_{GLT,end}')
% xlabel('v_{GLT,20}')
% ylabel('v_{GLT,end}')

xlabel('v_{HXT} at time = 20 s (mM s^{-1})')
ylabel('v_{HXT} at time = 400 s (mM s^{-1})')

% ylabel('v_{GLT,20}')
% xlabel('v_{GLT,end}')
% % % % ylabel('vHXT at 20 s (mM s^{-1})')
% % % % xlabel('vHXT at 0 or 400 s (mM s^{-1})')
% % % % ylabel('GLCec at 0 or 400 s (mM)')
% % % xlabel('error GLC')
% % % ylabel('error GLC')

xlim([0 2])%50])
ylim([0 0.4])
xticks([0 0.5 1 1.5 2])
yticks([0 0.1 0.2 0.3 0.4])

% xlim([0 0.4])%50])
% ylim([0 1.5])

% specific data point
hold on
% plot(dataset.FF01.fluxes{1}(5), dataset.FF01.fluxes{1}(2), ...
%     'ko','MarkerFaceColor','b')
plot(dataset.FF01.fluxes{1}(2), dataset.FF01.fluxes{1}(5), ...
    'ko','MarkerFaceColor','r')
% plot(dataset.FF01.fluxes{1}(5), dataset.FF01.metabolites.ECglucose.conc(1), ...
%     'ko','MarkerFaceColor','b')
% legend('MPSA range', 'MPSA data points','experimental')

text(0.1, 0.36, 'A', 'FontSize', 20)

box on

% % 
% subplot(3,2,5)
% hold on
% for i = 1:lenS%100%1000
%     temp_x = simRes1{i}.T_FF01;
%     temp_y = simRes1{i}.Y_FF01(:,36);
%     if(isreal(temp_y) == 0)
%         continue
%     end
%     plot(temp_x, temp_y, ...
%         'k-', 'Linewidth', 0.5)
% end
% plot(simRes_initial{1}.T_FF01, simRes_initial{1}.Y_FF01(:,36), ...
%     'b-', 'Linewidth', 1.5)
% % 
% s = scatter(dataset.FF01.metabolites.ECglucose.time,...
%     dataset.FF01.metabolites.ECglucose.conc);
% s.LineWidth = 0.6;
% s.SizeData = 20;
% s.MarkerEdgeColor = 'k';
% s.MarkerFaceColor = 'r';
% ylabel('GLCec concentration (mM)')
% xlabel('time (s)')
% % ylim([0 0.6])
% % xlim([0 400])
% box on

% 
% % % % subplot(3,2,6)
    axes('Position',[.525 .175 .35 .35])
    box on

hold on
% plots in grey
n_temp = 200; %lenS
for i = 1:n_temp%lenS%100%1000
    temp_x = simRes1{i}.T_FF01;
    temp_y = simRes1{i}.V_FF01(:,1);
%     if(isreal(temp_y) == 0)
    if((isreal(temp_y) == 0)||(max(temp_y)) > 2)
        continue
    end
    % 
    if mod(i,21) == 0
%         plot(temp_x, temp_y, ...
%             'k-', 'Linewidth', 0.5, ...
%             'Color', [.15 .15 .15])
    else
        plot(temp_x, temp_y, ...
            '-', 'Linewidth', 0.5, ...
            'Color', [.85 .85 .85])
    end
end
% plots in black
for i = 1:n_temp%lenS%100%1000
    temp_x = simRes1{i}.T_FF01;
    temp_y = simRes1{i}.V_FF01(:,1);
%     if(isreal(temp_y) == 0)
    if((isreal(temp_y) == 0)||(max(temp_y)) > 2)
        continue
    end
    % 
    if mod(i,21) == 0
        plot(temp_x, temp_y, ...
            'k-', 'Linewidth', 0.5, ...
            'Color', [.15 .15 .15])
    else
%         plot(temp_x, temp_y, ...
%             '-', 'Linewidth', 0.5, ...
%             'Color', [.85 .85 .85])
    end
end
% nice fit
plot(simRes_initial{1}.T_FF01, simRes_initial{1}.V_FF01(:,1), ...
    'b-', 'Linewidth', 1.5)
% experimental data
s = scatter(dataset.FF01.fluxes_times,...
    dataset.FF01.fluxes{1});
s.LineWidth = 0.6;
s.SizeData = 20;
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'r';
% ylabel('vHXT (mM s^{-1})')
% xlabel('time (s)')
xlab = xlabel('v_{HXT} vs time (mM s^{-1}, s)');
xlab.Position = xlab.Position + [0 2 0];
% ylim([0 3])
% xlim([0 400])

text(325, 1.25, 'B', 'FontSize', 20)

box on

% % 
set(10, 'Position', [100 100 600 500])
set(10,'Color','w');


% % % % % %% Saving outcome
% % % % % set(gcf,'Units','inches');
% % % % % screenposition = get(gcf,'Position');
% % % % % set(gcf,...
% % % % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % % % %     'PaperSize',[screenposition(3:4)]);
% % % % % print -dpdf -painters Y3M2_glt_csmin
% % % % % % print -depsc -painters figure_2
% % % % % % print('figure_2','-depsc','-tiff')
% % % % % print -dpng -painters Y3M2_glt_csmin


%%
% % % % % clf(11)
figure(11)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
subplot(3,1,1)
% scatter data
hold on
% 
% parsGLT = [35 36 38];
s = scatter(vGLT_20,vGLT_end,20, xMPSA(:,35),'Filled');
% 
s.LineWidth = 0.6;
s.SizeData = 15;
% 
xlabel('v_{HXT} at time = 20 s (mM s^{-1})')
ylabel('v_{HXT} at time = 400 s (mM s^{-1})')
% 
xlim([0 1.5])%50])
ylim([0 0.4])
xticks([0 0.5 1 1.5])
yticks([0 0.1 0.2 0.3 0.4])
% specific data point
hold on
plot(dataset.FF01.fluxes{1}(2), dataset.FF01.fluxes{1}(5), ...
    'ko','MarkerFaceColor','r')

text(0.1, 0.36, 'Km.GLTic', 'FontSize', 14)

box on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
subplot(3,1,2)
% scatter data
hold on
% 
% parsGLT = [35 36 38];
s = scatter(vGLT_20,vGLT_end,20, xMPSA(:,36),'Filled');
% 
s.LineWidth = 0.6;
s.SizeData = 15;
% 
xlabel('v_{HXT} at time = 20 s (mM s^{-1})')
ylabel('v_{HXT} at time = 400 s (mM s^{-1})')
% 
xlim([0 1.5])%50])
ylim([0 0.4])
xticks([0 0.5 1 1.5])
yticks([0 0.1 0.2 0.3 0.4])
% specific data point
hold on
plot(dataset.FF01.fluxes{1}(2), dataset.FF01.fluxes{1}(5), ...
    'ko','MarkerFaceColor','r')

text(0.1, 0.36, 'Vmax', 'FontSize', 14)

box on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
subplot(3,1,3)
% scatter data
hold on
% 
% parsGLT = [35 36 38];
s = scatter(vGLT_20,vGLT_end,20, xMPSA(:,38),'Filled');
% 
s.LineWidth = 0.6;
s.SizeData = 15;
% 
xlabel('v_{HXT} at time = 20 s (mM s^{-1})')
ylabel('v_{HXT} at time = 400 s (mM s^{-1})')
% 
xlim([0 1.5])%50])
ylim([0 0.4])
xticks([0 0.5 1 1.5])
yticks([0 0.1 0.2 0.3 0.4])
% specific data point
hold on
plot(dataset.FF01.fluxes{1}(2), dataset.FF01.fluxes{1}(5), ...
    'ko','MarkerFaceColor','r')

text(0.1, 0.36, 'Km.GLTec', 'FontSize', 14)

box on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % 
set(11, 'Position', [100 100 300 900])
set(11,'Color','w');


% % % % % %% Saving outcome
% % % % % set(gcf,'Units','inches');
% % % % % screenposition = get(gcf,'Position');
% % % % % set(gcf,...
% % % % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % % % %     'PaperSize',[screenposition(3:4)]);
% % % % % print -dpdf -painters Y3M2_glt_csmin_sup
% % % % % % print -depsc -painters figure_3
% % % % % % print('figure_2','-depsc','-tiff')
% % % % % print -dpng -painters Y3M2_glt_csmin_sup

