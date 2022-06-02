% % 
% Figura A in 'embed simulations'

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
setup.clamp_GLCec = 0;
% setup.clamp_GLCec = 1;
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


%% (2021 10 01: adding the sinkPYR decrease)
setup.decrease_sinkPYR = 1;
setup.ratio_decrease_sinkPYR = 0.1;


%% Starting simulation
setup.experiment = 1;
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% plotMode = 2; % multiple simulations
% % 
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;
% % % 
plotMode = 2;
referencePlotSimulations_enrichment
plotMode = 0;

% 
simRes_initial = selResCell;
clear selResCell
% selResCell = simRes_initial;

%% 
% Incoming reactions
v_in_time = simRes_initial{1}.T_FF01;
v_in_GLT = simRes_initial{1}.V_FF01(:,1); % + v_GLT 
v_in_NTH1 = 2*simRes_initial{1}.V_FF01(:,20); % + 2.*v_NTH1 
v_in_ATH1v = 2*simRes_initial{1}.V_FF01(:,46); % + 2.*v_ATH1vac 
v_in_GLYCdeg = simRes_initial{1}.V_FF01(:,51); % + v_glycDeg
% Labeling percentages
perc_lab_GLCec = simRes_initial{1}.Y_FF01(:,65) ./ simRes_initial{1}.Y_FF01(:,36); % GLCec
perc_lab_TREcyt = simRes_initial{1}.Y_FF01(:,58) ./ simRes_initial{1}.Y_FF01(:,25); % TREcyt
perc_lab_TREvac = simRes_initial{1}.Y_FF01(:,67) ./ simRes_initial{1}.Y_FF01(:,38); % TREvac
perc_lab_GLYCcyt = simRes_initial{1}.Y_FF01(:,68) ./ simRes_initial{1}.Y_FF01(:,40); % GLYCcyt
% % Incoming fluxes labelled
v_in_GLT_lab = v_in_GLT .* perc_lab_GLCec;
v_in_NTH1_lab = v_in_NTH1 .* perc_lab_TREcyt;
v_in_ATH1v_lab = v_in_ATH1v .* perc_lab_TREvac;
v_in_GLYCdeg_lab = v_in_GLYCdeg .* perc_lab_GLYCcyt;
% % Incoming fluxes unlabelled
v_in_GLT_unlab = v_in_GLT - v_in_GLT_lab;
v_in_NTH1_unlab = v_in_NTH1 - v_in_NTH1_lab;
v_in_ATH1v_unlab = v_in_ATH1v - v_in_ATH1v_lab;
v_in_GLYCdeg_unlab = v_in_GLYCdeg - v_in_GLYCdeg_lab;

% % Plot all first, each category per row
% clf(21)
% figure(21)
if exist('fh_21','var')
    clf(21)
else
    fh_21 = figure(21);
end
% 
subplot(4,4,1), plot(v_in_time,v_in_GLT,'k-'), title('v, GLT')
subplot(4,4,2), plot(v_in_time,v_in_NTH1,'k-'), title('v, NTH1')
subplot(4,4,3), plot(v_in_time,v_in_ATH1v,'k-'), title('v, ATH1v')
subplot(4,4,4), plot(v_in_time,v_in_GLYCdeg,'k-'), title('v, GLYCdeg')
% 
subplot(4,4,5), plot(v_in_time,perc_lab_GLCec,'k-'), title('perc, GLCec')
subplot(4,4,6), plot(v_in_time,perc_lab_TREcyt,'k-'), title('perc, TREcyt')
subplot(4,4,7), plot(v_in_time,perc_lab_TREvac,'k-'), title('perc, TREvac')
subplot(4,4,8), plot(v_in_time,perc_lab_GLYCcyt,'k-'), title('perc, GLYCcyt')
% 
subplot(4,4,9), plot(v_in_time,v_in_GLT_lab,'k-'), title('lab, GLT')
subplot(4,4,10), plot(v_in_time,v_in_NTH1_lab,'k-'), title('lab, NTH1')
subplot(4,4,11), plot(v_in_time,v_in_ATH1v_lab,'k-'), title('lab, ATH1v')
subplot(4,4,12), plot(v_in_time,v_in_GLYCdeg_lab,'k-'), title('lab, GLYCdeg')
% 
subplot(4,4,13), plot(v_in_time,v_in_GLT_unlab,'k-'), title('unlab, GLT')
subplot(4,4,14), plot(v_in_time,v_in_NTH1_unlab,'k-'), title('unlab, NTH1')
subplot(4,4,15), plot(v_in_time,v_in_ATH1v_unlab,'k-'), title('unlab, ATH1v')
subplot(4,4,16), plot(v_in_time,v_in_GLYCdeg_unlab,'k-'), title('unlab, GLYCdeg')
% 

%% making the area plot
% https://colorbrewer2.org/#type=sequential&scheme=PuBu&n=5
% % 
% col1 = [255, 87, 51]/255;
% col2 = [239, 144, 124]/255;
% col3 = [213, 105, 81]/255;
% col4 = [183, 71, 47]/255;
% col5 = [47, 69, 183]/255;
% col6 = [110, 129, 233]/255;
% col7 = [102, 113, 174]/255;
% col8 = [44, 64, 170]/255;
% (for 5, select 4 top)
col1 = [165,15,21]/255;
col2 = [222,45,38]/255;
col3 = [251,106,74]/255;
col4 = [252,174,145]/255;
% 
col8 = [189,215,231]/255;
col7 = [107,174,214]/255;
col6 = [49,130,189]/255;
col5 = [8,81,156]/255;
%
colArr = [col1; col2; col3; col4; ...
    col5; col6; col7; col8];
% 
if exist('enrichmentArea','var')
    clf(22)
end
% amounts vs time
enrichmentArea = figure(22);
enrichmentArea.Position = [1954 374 237 455];
hold on
sp1 = subplot(2,1,1);
temp_area = area(v_in_time,[v_in_GLT_lab, v_in_NTH1_lab, v_in_ATH1v_lab, v_in_GLYCdeg_lab, ...
                v_in_GLT_unlab, v_in_NTH1_unlab, v_in_ATH1v_unlab, v_in_GLYCdeg_unlab]);
% 
for i = 1:8
    temp_area(i).EdgeColor = 'none';
    temp_area(i).FaceColor = colArr(i,:);
end

% percentages over total
v_in_total = v_in_GLT_lab + v_in_NTH1_lab + v_in_ATH1v_lab + v_in_GLYCdeg_lab + ...
                v_in_GLT_unlab + v_in_NTH1_unlab + v_in_ATH1v_unlab + v_in_GLYCdeg_unlab;
% 
ylim([0 0.6])
yticks([0 0.2 0.4 0.6])
xlim([0 400])
xticks([0 100 200 300 400])
%
xlabel('Time (s)')
ylabel('Inflow GLC (mM s^{-1})')


% 
sp2 = subplot(2,1,2);
hold on
temp_area_perc = area(v_in_time,[v_in_GLT_lab./v_in_total, v_in_NTH1_lab./v_in_total, v_in_ATH1v_lab./v_in_total, v_in_GLYCdeg_lab./v_in_total, ...
                v_in_GLT_unlab./v_in_total, v_in_NTH1_unlab./v_in_total, v_in_ATH1v_unlab./v_in_total, v_in_GLYCdeg_unlab./v_in_total]);
% % 
% legend('GLT_{lab}', 'NTH1_{lab}', 'ATH1v_{lab}', 'GLYCdeg_{lab}', ...
%     'GLT_{unlab}', 'NTH1_{unlab}', 'ATH1v_{unlab}', 'GLYCdeg_{unlab}')
% 
for i = 1:8
    temp_area_perc(i).EdgeColor = 'none';
    temp_area_perc(i).FaceColor = colArr(i,:);
end
v_temp = v_in_GLT_lab./v_in_total + v_in_NTH1_lab./v_in_total + v_in_ATH1v_lab./v_in_total + v_in_GLYCdeg_lab./v_in_total;

% box on
plot(v_in_time, v_temp, 'k-', 'LineWidth', 0.5)
% 
ylim([0 1])
yticks([0 0.25 0.5 0.75 1])
xlim([0 400])
xticks([0 100 200 300 400])
%
xlabel('Time (s)')
ylabel('Inflow GLC (%)')

% last edits: legend
temp_leg = legend('GLT_{lab}', 'NTH1_{lab}', 'ATH1v_{lab}', 'GLYCdeg_{lab}', ...
    'GLT_{unlab}', 'NTH1_{unlab}', 'ATH1v_{unlab}', 'GLYCdeg_{unlab}');
% % horizontal tab
% temp_leg = legend('GLT_{lab}', 'NTH1_{lab}', 'ATH1v_{lab}', 'GLYCdeg_{lab}', ...
%     'GLT_{unlab}', 'NTH1_{unlab}', 'ATH1v_{unlab}', 'GLYCdeg_{unlab}',...
%     'Orientation','Horizontal');

% % % % % %% Saving
% % % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F7CD');
% % % % % %% only for the legend
% % % % % set(enrichmentArea,'Units','Inches');
% % % % % pos = get(enrichmentArea,'Position');
% % % % % set(enrichmentArea,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % % print(enrichmentArea,'E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F7CDleg','-dpdf','-r0')


