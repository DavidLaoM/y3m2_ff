% % JUST RUN BELOW FOR NOW


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
x2 = x;
% % % % % % % % % % % %% adjustment glk regularization
% % % % % % % % % % parsGLK = [28, 29, 30, 31, 32, 33, 34]; 
% % % % % % % % % % load('FF_pE16a2_glk_wA10.mat');
% % % % % % % % % % x(parsGLK) = xres;


%% Directly implementing the Csmin idea in out model
% general clamping and setup options
NumberCycles = 20; % 5 %20;
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
% 
setup.experiment = 1;


%% (2021 10 01: adding the sinkPYR decrease)
setup.decrease_sinkPYR = 1;
setup.ratio_decrease_sinkPYR = 0.1;


% %% Simulation
% % 
% [T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % 
% % T_FF01_1 = T_FF01_2;
% % Y_FF01_1 = Y_FF01_2;
% % V_FF01_1 = V_FF01_2;
% % 
% selResCell{1}.T_FF01 = T_FF01_1;
% selResCell{1}.Y_FF01 = Y_FF01_1;
% selResCell{1}.V_FF01 = V_FF01_1;


% %% Visualization
% plotMode = 2;
% referencePlotSimulations_enrichment
% % clear selResCell
% plotMode = 0;


%% Check simulation with the latest regularization
% % % % % % % % % % % 
% % % % % % % % % % load('parSet_105b.mat','x105b');
% % % % % % % % % % x2 = x;
% % % % % % % % % %     parsPFK = 43:56;
% % % % % % % % % % %     parsALD = 11:15;
% % % % % % % % % %     parsPYK = 71:77;
% % % % % % % % % % x2(parsPFK) = x105b(parsPFK);
% % % % % % % % % % % x2(parsALD) = x105b(parsALD);
% % % % % % % % % % x2(parsPYK) = x105b(parsPYK);
% % % % % % % % % % x2(152) = -2.4128;
% % % % % % % % % % % 
% % % % % % % % % % x2 = x;
% 
[T_FF01_2,Y_FF01_2,V_FF01_2] = simulate_FF_enrichment(x2,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % 
% % 
% selResCell{2}.T_FF01 = T_FF01_2;
% selResCell{2}.Y_FF01 = Y_FF01_2;
% selResCell{2}.V_FF01 = V_FF01_2;
% 
selResCell{1}.T_FF01 = T_FF01_2;
selResCell{1}.Y_FF01 = Y_FF01_2;
selResCell{1}.V_FF01 = V_FF01_2;


%% Visualization
plotMode = 2;
referencePlotSimulations_enrichment
% clear selResCell
plotMode = 0;


%% plot in here
% disp('plot in here')
array_time = selResCell{1}.T_FF01';
array_blank = ones(size(array_time));
interp_time = 0:1:400;
interp_array_blank = ones(size(interp_time));

% % line UG
array_GLCext = (selResCell{1}.Y_FF01(:,36)' - min(selResCell{1}.Y_FF01(:,36)));
array_GLCext = array_GLCext./max(array_GLCext);
array_GLCext_int = interp1(array_time,array_GLCext,interp_time,'pchip');
%
array_GLC = (selResCell{1}.Y_FF01(:,6)' - min(selResCell{1}.Y_FF01(:,6)));
array_GLC = array_GLC./max(array_GLC);
array_GLC_int = interp1(array_time,array_GLC,interp_time,'pchip');
%
array_G6P = (selResCell{1}.Y_FF01(:,5)' - min(selResCell{1}.Y_FF01(:,5)));
array_G6P = array_G6P./max(array_G6P);
array_G6P_int = interp1(array_time,array_G6P,interp_time,'pchip');
%
array_F6P = (selResCell{1}.Y_FF01(:,4)' - min(selResCell{1}.Y_FF01(:,4)));
array_F6P = array_F6P./max(array_F6P);
array_F6P_int = interp1(array_time,array_F6P,interp_time,'pchip');
%
array_FBP = (selResCell{1}.Y_FF01(:,3)' - min(selResCell{1}.Y_FF01(:,3)));
array_FBP = array_FBP./max(array_FBP);
array_FBP_int = interp1(array_time,array_FBP,interp_time,'pchip');

% % line LG
array_GAP = (selResCell{1}.Y_FF01(:,14)' - min(selResCell{1}.Y_FF01(:,14)));
array_GAP = array_GAP./max(array_GAP);
array_GAP_int = interp1(array_time,array_GAP,interp_time,'pchip');
%
array_BPG = (selResCell{1}.Y_FF01(:,2)' - min(selResCell{1}.Y_FF01(:,2)));
array_BPG = array_BPG./max(array_BPG);
array_BPG_int = interp1(array_time,array_BPG,interp_time,'pchip');
%
array_3PG = (selResCell{1}.Y_FF01(:,11)' - min(selResCell{1}.Y_FF01(:,11)));
array_3PG = array_3PG./max(array_3PG);
array_3PG_int = interp1(array_time,array_3PG,interp_time,'pchip');
%
array_2PG = (selResCell{1}.Y_FF01(:,10)' - min(selResCell{1}.Y_FF01(:,10)));
array_2PG = array_2PG./max(array_2PG);
array_2PG_int = interp1(array_time,array_2PG,interp_time,'pchip');
%
array_PEP = (selResCell{1}.Y_FF01(:,12)' - min(selResCell{1}.Y_FF01(:,12)));
array_PEP = array_PEP./max(array_PEP);
array_PEP_int = interp1(array_time,array_PEP,interp_time,'pchip');
%
array_PYR = (selResCell{1}.Y_FF01(:,13)' - min(selResCell{1}.Y_FF01(:,13)));
array_PYR = array_PYR./max(array_PYR);
array_PYR_int = interp1(array_time,array_PYR,interp_time,'pchip');

% % line storage
array_G1P = (selResCell{1}.Y_FF01(:,21)' - min(selResCell{1}.Y_FF01(:,21)));
array_G1P = array_G1P./max(array_G1P);
array_G1P_int = interp1(array_time,array_G1P,interp_time,'pchip');
%
array_UDPG = (selResCell{1}.Y_FF01(:,24)' - min(selResCell{1}.Y_FF01(:,24)));
array_UDPG = array_UDPG./max(array_UDPG);
array_UDPG_int = interp1(array_time,array_UDPG,interp_time,'pchip');
%
array_T6P = (selResCell{1}.Y_FF01(:,26)' - min(selResCell{1}.Y_FF01(:,26)));
array_T6P = array_T6P./max(array_T6P);
array_T6P_int = interp1(array_time,array_T6P,interp_time,'pchip');
%
array_TREic = (selResCell{1}.Y_FF01(:,38)' - min(selResCell{1}.Y_FF01(:,38)));
array_TREic = array_TREic./max(array_TREic);
array_TREic_int = interp1(array_time,array_TREic,interp_time,'pchip');

% % line glycerol
array_DHAP = (selResCell{1}.Y_FF01(:,17)' - min(selResCell{1}.Y_FF01(:,17)));
array_DHAP = array_DHAP./max(array_DHAP);
array_DHAP_int = interp1(array_time,array_DHAP,interp_time,'pchip');
%
array_G3P = (selResCell{1}.Y_FF01(:,18)' - min(selResCell{1}.Y_FF01(:,18)));
array_G3P = array_G3P./max(array_G3P);
array_G3P_int = interp1(array_time,array_G3P,interp_time,'pchip');
%
array_GLYC = (selResCell{1}.Y_FF01(:,19)' - min(selResCell{1}.Y_FF01(:,19)));
array_GLYC = array_GLYC./max(array_GLYC);
array_GLYC_int = interp1(array_time,array_GLYC,interp_time,'pchip');


% yticks([0 50 100])
% yticklabels({'y = 0','y = 50','y = 100'})
yticks_vals = flip(linspace(0.5,16.5,17));
yticklabels_vals = {'GLCext','GLC','G6P','F6P','FBP',...
    'GAP','BPG','3PG','2PG','PEP','PYR',...
    'G1P','T6P','UDPG',...%'TRE',...
    'DHAP','G3P','GLYC'};


%%
colArr2 = zeros(401,3); colArr2(1:100,1) = 0.5;
% 
if exist('fh121','var')
    clf(121)
end
fh121 = figure(121);
% 
fill_area = 0:0.01:1;
for j = 1:length(yticks_vals)
    for i = fill_area % it's in the y value where it is plotted, then that strange range
        if j == 1, scatter(interp_time,16+interp_array_blank*i,1,array_GLCext_int,'filled'),
        elseif j == 2, scatter(interp_time,15+interp_array_blank*i,1,array_GLC_int,'filled'),
        elseif j == 3, scatter(interp_time,14+interp_array_blank*i,1,array_G6P_int,'filled'),
        elseif j == 4, scatter(interp_time,13+interp_array_blank*i,1,array_F6P_int,'filled'),
        elseif j == 5, scatter(interp_time,12+interp_array_blank*i,1,array_FBP_int,'filled'),
        elseif j == 6, scatter(interp_time,11+interp_array_blank*i,1,array_GAP_int,'filled'),
        elseif j == 7, scatter(interp_time,10+interp_array_blank*i,1,array_BPG_int,'filled'),
        elseif j == 8, scatter(interp_time,9+interp_array_blank*i,1,array_3PG_int,'filled'),
        elseif j == 9, scatter(interp_time,8+interp_array_blank*i,1,array_2PG_int,'filled'),
        elseif j == 10, scatter(interp_time,7+interp_array_blank*i,1,array_PEP_int,'filled'),
        elseif j == 11, scatter(interp_time,6+interp_array_blank*i,1,array_PYR_int,'filled'),
        elseif j == 12, scatter(interp_time,5+interp_array_blank*i,1,array_G1P_int,'filled'),
        elseif j == 13, scatter(interp_time,4+interp_array_blank*i,1,array_T6P_int,'filled'),
        elseif j == 14, scatter(interp_time,3+interp_array_blank*i,1,array_UDPG_int,'filled'),
%         elseif j == 15, scatter(interp_time,14+interp_array_blank*i,1,array_TREic_int,'filled'),
        elseif j == 15, scatter(interp_time,2+interp_array_blank*i,1,array_DHAP_int,'filled'),
        elseif j == 16, scatter(interp_time,1+interp_array_blank*i,1,array_G3P_int,'filled'),
        elseif j == 17, scatter(interp_time,interp_array_blank*i,1,array_GLYC_int,'filled'),
        end
        colormap(flipud(gray))
        hold on
    end
end
% 
line([0 0],[0 17],'color','k','LineWidth',1)
line([400 400],[0 17],'color','k','LineWidth',1)
% 
for i = 0:1:length(yticks_vals)
    line([0 400],[i i],'color','k','LineWidth',1)
%     line([0 400],[0 0],'color','k','LineWidth',1)
%     line([0 400],[1 1],'color','k','LineWidth',1)
%     line([0 400],[2 2],'color','k','LineWidth',1)
%     line([0 400],[3 3],'color','k','LineWidth',1)
end
% %%
yticks(flip(yticks_vals))
yticklabels(flip(yticklabels_vals))
xlim([0 400])
ylim([0 17])
xlab = xlabel('Time (s)');
ylab = ylabel('Metabolite concentration');
box on

% %% extra layout
% safecopy
fh121.InnerPosition = [2711 474 560 420];
fh121.OuterPosition = [2703 466 576 513];
fh121.Position = [2711 474 560 420];
% additions
% fh121.OuterPosition = [2703 466 576*0.9 513*0.9];
fh121.Position = [2711 474 560*1.1 420*1.1];
% fh121.InnerPosition = [2711 474 560*0.9 420*0.9];
fh121_children = get(fh121, 'Children');
fh121_children.Position = [0.1324+0.04 0.1100+0.02 0.7726-0.08 0.8150-0.04];
% plot([0 400],[0 17],'r-')
% hold off
% plot([0 20], [0 17.5], 'k-', 'Linewidth', 5)
% % % % temp = text(2,8,'A Simple Plot','Color','red','FontSize',14)


% %% top labeling
dim = [.177 .93 .025 .001];
str = '';
temp = annotation('textbox',dim,'String',str,'FitBoxToText','on', 'LineWidth', 5);
temp2 = annotation('textbox',[.177+.025 .93 0.6681 .001],'String',str,...
    'FitBoxToText','on', 'LineWidth', 1,'LineStyle','--');

% %%
% tempFe = text(10,18.5,'Feast','Color','black','FontSize',10,'HorizontalAlignment','center');
% tempFe = text(10,18.5,'Feast','Color','black','FontSize',10,'HorizontalAlignment','center');
% tempFa = text(190+20,18.5,'Famine','Color','black','FontSize',10,'HorizontalAlignment','center');
% tempFa = text(190+20,18.5,'Famine','Color','black','FontSize',10,'HorizontalAlignment','center');
tempFe = text(10,18.5,'Feeding','Color','black','FontSize',10,'HorizontalAlignment','center');
tempFe = text(10,18.5,'Feeding','Color','black','FontSize',10,'HorizontalAlignment','center');
tempFa = text(190+20,18.5,'No Feeding','Color','black','FontSize',10,'HorizontalAlignment','center');
tempFa = text(190+20,18.5,'No Feeding','Color','black','FontSize',10,'HorizontalAlignment','center');


% %% left labeling
% tempGLYC = annotation('textbox',[.177+.025+0.68 .775 .001 0.12],'String',str,'FitBoxToText','on', 'LineWidth', 1);
tempGLYC = annotation('textbox',[.177+.025+0.68 .775-.135-.230-0.2725 .001 0.12],'String',str,'FitBoxToText','on', 'LineWidth', 1);
% %%
% delete(tempGLYC)
% %%
% tempSTOR = annotation('textbox',[.177+.025+0.68 .775-.135 .001 0.12],'String',str,'FitBoxToText','on', 'LineWidth', 1);
tempSTOR = annotation('textbox',[.177+.025+0.68 .775-.230-0.2725 .001 0.12],'String',str,'FitBoxToText','on', 'LineWidth', 1);
% %%
% delete(tempSTOR)
% %%
% tempLG = annotation('textbox',[.177+.025+0.68 .775-.135-.230 .001 0.12+0.0875],'String',str,'FitBoxToText','on', 'LineWidth', 1);
tempLG = annotation('textbox',[.177+.025+0.68 .775-.230-0.2725+.140 .001 0.12+0.0875],'String',str,'FitBoxToText','on', 'LineWidth', 1);
% %%
% delete(tempLG)
% %%
% tempUG = annotation('textbox',[.177+.025+0.68 .775-.135-.230-0.2725 .001 0.12+0.0875+0.04],'String',str,'FitBoxToText','on', 'LineWidth', 1);
tempUG = annotation('textbox',[.177+.025+0.68 .775-0.2725+.140 .001 0.12+0.0875+0.04],'String',str,'FitBoxToText','on', 'LineWidth', 1);
% %%
% delete(tempUG)

% %%
textGLYC = text(425,1.5,'Glycerol','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);
textGLYC = text(425,1.5,'Glycerol','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);
% %%
textTRE = text(450,1.5+3,'Trehalose','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);
textTRE = text(450,1.5+3,'Trehalose','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);
% %%
textLG = text(425,1.5+3+4,'Lower Glycolysis','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);
textLG = text(425,1.5+3+4,'Lower Glycolysis','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);
% %%
textUG = text(450,1.5+3+4+5.5,'Upper Glycolysis','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);
textUG = text(450,1.5+3+4+5.5,'Upper Glycolysis','Color','black','FontSize',10,'HorizontalAlignment','center','Rotation',90);

% %%
% delete(textUG)

% %% color babckground
set(fh121,'color','w')
%
a = annotation('textbox',[.05 .7 .3 .3],'String','C',...
    'FitBoxToText','on','FontSize',20,'EdgeColor','none');


%%

% % % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F4Cbackbone');

