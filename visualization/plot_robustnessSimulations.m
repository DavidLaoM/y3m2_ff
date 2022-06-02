%% (to be located outside in the end) Initial setup:
robust_setup.simulation_reference = selResIntial;
robust_setup.simulation_robustness = simRobust;
%     robust_setup.simulation_robustness = selResIntial;
%     for i = 1:39
%         robust_setup.simulation_robustness = [robust_setup.simulation_robustness, selResIntial];
%     end
robust_setup.forced_density = 0;
% robust_setup.simulation_robustness_idxs = 1:10; 
robust_setup.simulation_robustness_idxs = 1:1000; %1:100 %1:40
% robust_setup.simulation_robustness_idxs = 1:10000; %1:100 %1:40
% robust_setup.simulation_robustness_idxs = 1:3000; %1:100 %1:40
% % % % robust_setup.simulation_robustness_idxs = 1:100; %1:40
robust_setup.interpSpace = [0:400]';
% robust_setup.dpVal = 256;
robust_setup.dpVal = 400;
% robust_setup.dpVal_yaxis = 40;
robust_setup.dpVal_yaxis = 75;
% robust_setup.dpVal_yaxis = 400;
robust_setup.maxRange_time = 400;
robust_setup.LineWidth = 2;
robust_setup.MarkerSize = 30;
robust_setup.contourLevels = 500; %32
robust_setup.Ylim_foldIncrease = 3;
robust_setup.shading_increase_ratio = 0; % 0.1

% metabolite-specific
% % % % robust_setup.metabolite_id = 15;
% % % % robust_setup.maxRange_concentration = 5;
robust_setup.nMets = 38;
robust_setup.metabolite_id_range = 1:38; %1:15; % 1:38
robust_setup.maxRange_concentration_range = zeros(1,robust_setup.nMets);
% loop to get the ymax from experimental data and, if not, from the
% simulation
f1_h = figure(1);
f1_c_h = get(f1_h,'children');
% 
for i = 1:length(robust_setup.maxRange_concentration_range)
    % option based on simulated axes elsewhere
    robust_setup.maxRange_concentration_range(i) = f1_c_h(end+1-i).YLim(2);
%     Another option
%     if isfield(ExpData.metabolites{i}, 'conc')
%         robust_setup.maxRange_concentration_range(i) = max(ExpData.metabolites{i}.conc)*robust_setup.Ylim_foldIncrease;
%     else
%         robust_setup.maxRange_concentration_range(i) = max(selResIntial{1}.Y_FF01(:,i))*robust_setup.Ylim_foldIncrease;
%     end
end
% 
robust_setup.maxRange_concentration_range(1) = 0.5;
% robust_setup.maxRange_concentration_range(2) = 0.01;
% robust_setup.maxRange_concentration_range(3) = 2;
% robust_setup.maxRange_concentration_range(4) = 2;
% robust_setup.maxRange_concentration_range(5) = 10;
% robust_setup.maxRange_concentration_range(6) = 0.4;
% 
% robust_setup.maxRange_concentration_range(7) = 2;
% robust_setup.maxRange_concentration_range(8) = 0.5;
robust_setup.maxRange_concentration_range(9) = 8;
% robust_setup.maxRange_concentration_range(10) = 1;
% robust_setup.maxRange_concentration_range(11) = 10;
% robust_setup.maxRange_concentration_range(12) = 4;
% 
% robust_setup.maxRange_concentration_range(13) = 4;
% robust_setup.maxRange_concentration_range(14) = 0.04;
robust_setup.maxRange_concentration_range(15) = 2.5;
% robust_setup.maxRange_concentration_range(16) = 1;
% robust_setup.maxRange_concentration_range(17) = 1;
% robust_setup.maxRange_concentration_range(18) = 0.1;
% 
robust_setup.maxRange_concentration_range(19) = 0.5;
% robust_setup.maxRange_concentration_range(20) = 20;
% robust_setup.maxRange_concentration_range(21) = 0.4;
% robust_setup.maxRange_concentration_range(22) = 1;
% robust_setup.maxRange_concentration_range(23) = 0.4;
% robust_setup.maxRange_concentration_range(24) = 4;
% 
% robust_setup.maxRange_concentration_range(25) = 10;
% robust_setup.maxRange_concentration_range(26) = 2;
robust_setup.maxRange_concentration_range(27) = 10.1;
% robust_setup.maxRange_concentration_range(28) = 0.05;
% robust_setup.maxRange_concentration_range(29) = 0.05;
% robust_setup.maxRange_concentration_range(30) = 1;
% 
robust_setup.maxRange_concentration_range(31) = 2.1;
robust_setup.maxRange_concentration_range(32) = 0.1;
robust_setup.maxRange_concentration_range(33) = 1;
robust_setup.maxRange_concentration_range(34) = 1;
robust_setup.maxRange_concentration_range(35) = 1;
robust_setup.maxRange_concentration_range(36) = 1;
% 
% robust_setup.maxRange_concentration_range(37) = 0.02;
robust_setup.maxRange_concentration_range(38) = 80;


% reaction rate-specific
robust_setup.nRates = 51;
robust_setup.rate_id_range = 1:51;
robust_setup.maxRange_rate_range = zeros(1,robust_setup.nRates);
% loop to get the ymax from experimental data and, if not, from the
% simulation
f2_h = figure(2);
f2_c_h = get(f2_h,'children');
% 
for i = 1:length(robust_setup.maxRange_rate_range)
    % option based on simulated axes elsewhere
    robust_setup.maxRange_rate_range(i) = f2_c_h(end+1-i).YLim(2);
end
% 
robust_setup.maxRange_rate_range(6) = 0.1;
% 
robust_setup.maxRange_rate_range(15) = 0.1;
% 
robust_setup.maxRange_rate_range(24) = 0.1;
% 
robust_setup.maxRange_rate_range(30) = 0.1;
robust_setup.maxRange_rate_range(31) = 0.1;
robust_setup.maxRange_rate_range(32) = 0.1;
% 
% 
robust_setup.maxRange_rate_range(42) = 0.1;
robust_setup.maxRange_rate_range(43) = 0.1;
robust_setup.maxRange_rate_range(44) = 0.1;

%% Figure for metabolites
if exist('fig_h','var')
    clear fig_h
end
fig_h = figure(1004);
for j = robust_setup.metabolite_id_range
    if j ~= 1
        toc
    end
    disp(j)
    tic
    % minimum setup each loop
    robust_setup.metabolite_id = j;
    robust_setup.maxRange_concentration = robust_setup.maxRange_concentration_range(j);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % %% arrange data
    % heatmap data
    if exist('tData','var')
        clear tData yData
    end
    tData = [];
    yData = [];
% % % %     for i = 1:robust_setup.simulation_robustness_idxs
    for i = robust_setup.simulation_robustness_idxs
        tData = [tData; robust_setup.interpSpace];
% % % %         yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.Y_FF01(:,15), robust_setup.interpSpace, 'pchip')];
        yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.Y_FF01(:,j), robust_setup.interpSpace, 'pchip')];
    end
    tData = [0; robust_setup.maxRange_time; tData];
    yData = [0; robust_setup.maxRange_concentration; yData];

    % simulation data
    sim_time = robust_setup.simulation_reference{1}.T_FF01 * robust_setup.dpVal / max(tData); %[1 2 10 100];
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,15) * robust_setup.dpVal / max(yData);
    sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,j) * robust_setup.dpVal / max(yData);

    % experimental data
    if isfield(ExpData.metabolites{j}, 'conc')
        exp_time = ExpData.metabolites{robust_setup.metabolite_id}.time * robust_setup.dpVal / max(tData);
        exp_conc = robust_setup.dpVal - ExpData.metabolites{robust_setup.metabolite_id}.conc * robust_setup.dpVal / max(yData);
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % %% Visualization (in loop)
    % % % % if exist('fig_h','var')
    % % % %     clear fig_h
    % % % % end
    % % % % fig_h = figure(1004);
    figure(fig_h);

    % make the heatmap plot
    sp_h = subplot(7,6,robust_setup.metabolite_id);
%     [f_h] = DataDensityPlot(tData, yData, robust_setup.contourLevels, fig_h, sp_h);
    [f_h, range90] = DataDensityPlot_1D(tData, yData, robust_setup.contourLevels, fig_h, sp_h, robust_setup);
    hold on

    % simulation data
% % % %     plot(sim_time, sim_conc,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/(max(sim_conc)-min(sim_conc))*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/max(sim_conc)*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
    plot(sim_time, sim_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)

%     % 
%     figure,
%     subplot(131), plot(range90(1,:),'.-')
%     subplot(132), plot(range90(2,:),'.-')
%     subplot(133), plot(range90(3,:),'.-')
    temp_range90min = robust_setup.dpVal - range90(2,:) * robust_setup.dpVal / max(yData);
    temp_range90max = robust_setup.dpVal - range90(3,:) * robust_setup.dpVal / max(yData);
    plot(range90(1,:), temp_range90min/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    plot(range90(1,:), temp_range90max/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,j) * robust_setup.dpVal / max(yData);
    
    
    % experimental data
    if isfield(ExpData.metabolites{j}, 'conc')
% % % %         scatter(exp_time, exp_conc, ...
        scatter(exp_time, exp_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
                    robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
                    'MarkerFaceColor', 'white', ...
                    'LineWidth', 1.5)
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % add text label
    metName_pre = erase(legenda.metabolites{j}, ", [mM]");
    metName = erase(metName_pre(1:end-7),",");
    if j == 14
        metName = 'GAP';
    elseif j == 19
        metName = 'GLYC';
    end
    text(sp_h.XLim(2)*0.7,sp_h.YLim(2)*0.25,metName,'FontSize',10)

end       


%% Figure for reaction rates
robust_setup.maxRange_concentration = 0; % for plots of rates, rathe than concentrations
if exist('fig_h2','var')
    clear fig_h2
end
% %%
fig_h2 = figure(1005);
for j = robust_setup.rate_id_range
    % minimum setup each loop
    robust_setup.rate_id = j;
    robust_setup.maxRange_rate = robust_setup.maxRange_rate_range(j);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % %% arrange data
    % heatmap data
    if exist('tData','var')
        clear tData yData
    end
    tData = [];
    yData = [];
% % % %     for i = 1:robust_setup.simulation_robustness_idxs
    for i = robust_setup.simulation_robustness_idxs
        tData = [tData; robust_setup.interpSpace];
% % % %         yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.Y_FF01(:,15), robust_setup.interpSpace, 'pchip')];
        yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.V_FF01(:,j), robust_setup.interpSpace, 'pchip')];
    end
    tData = [0; robust_setup.maxRange_time; tData];
    yData = [0; robust_setup.maxRange_rate; yData];

    % simulation data
    sim_time = robust_setup.simulation_reference{1}.T_FF01 * robust_setup.dpVal / max(tData); %[1 2 10 100];
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,15) * robust_setup.dpVal / max(yData);
    sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.V_FF01(:,j) * robust_setup.dpVal / max(yData);

    % experimental data
    if isfield(ExpData.fluxes{j}, 'rate')
        exp_time = ExpData.fluxes{robust_setup.rate_id}.time * robust_setup.dpVal / max(tData); exp_time(1) = 0.5;
        exp_flux = robust_setup.dpVal - ExpData.fluxes{robust_setup.rate_id}.rate * robust_setup.dpVal / max(yData);
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % %% Visualization (in loop)
    % % % % if exist('fig_h','var')
    % % % %     clear fig_h
    % % % % end
    % % % % fig_h = figure(1004);
    figure(fig_h2);

    % make the heatmap plot
    sp_h2 = subplot(7,8,robust_setup.rate_id);
%     [f_h] = DataDensityPlot(tData, yData, robust_setup.contourLevels, fig_h2, sp_h2);
    [f_h, range90] = DataDensityPlot_1D(tData, yData, robust_setup.contourLevels, fig_h2, sp_h2, robust_setup);
    hold on

    % simulation data
% % % %     plot(sim_time, sim_conc,'k-','LineWidth', robust_setup.LineWidth)
    plot(sim_time, sim_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)

    % 90% range
    temp_range90min = robust_setup.dpVal - range90(2,:) * robust_setup.dpVal / max(yData);
    temp_range90max = robust_setup.dpVal - range90(3,:) * robust_setup.dpVal / max(yData);
    plot(range90(1,:), temp_range90min/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    plot(range90(1,:), temp_range90max/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    
    % experimental data
    if isfield(ExpData.fluxes{j}, 'rate')
% % % %         scatter(exp_time, exp_conc, ...
        scatter(exp_time, exp_flux/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
                    robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
                    'MarkerFaceColor', 'white', ...
                    'LineWidth', 1.5)
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % add text label
    rateName_pre = erase(legenda.fluxes{j}, ", [mM s^{-1}]");
    rateName = erase(rateName_pre(1:end-7),",");
    text(sp_h.XLim(2)*0.7,sp_h.YLim(2)*0.25,rateName,'FontSize',10)
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
end       


% % % % %% safecopy 2
% % % % % Simulations with exprimental?tracing? data % % % % % % % % % % % % % % % % % % % 
% % % % 
% % % % 
% % % % figure%(2001)
% % % % for i = 1:38
% % % % 
% % % %     %FF01
% % % % 
% % % %     a1 = ExpData.metabolites{i};
% % % % %         t_sim = T_FF01;
% % % % %         c_sim = Y_FF01(:,i);
% % % %     if ~isempty(a1)
% % % %         tempFF01time = a1.time;
% % % %         tempFF01conc = a1.conc;
% % % % %             cip = tempFF01conc(tempFF01time<400);
% % % % %             tip = tempFF01time(tempFF01time<400);
% % % %         if isfield(a1, 'stdev')
% % % %             tempFF01std = a1.stdev;
% % % %         else 
% % % %             tempFF01std = zeros(size(tempFF01time)); 
% % % %         end
% % % % %             c_r2 = interp1(t_sim, c_sim, tip);
% % % % %             SSE = sum((cip - c_r2).^2);
% % % % %             SST = sum((cip - mean(cip)).^2);
% % % % %             R2 = 1 - SSE/SST;
% % % % %             R_matrix = corrcoef(c_r2, cip);
% % % % %             R = R_matrix(1,2);
% % % %     else
% % % %         tempFF01time = zeros(size(ExpData.metabolites{3}.time));
% % % %         tempFF01conc = zeros(size(ExpData.metabolites{3}.conc));
% % % %         tempFF01std = zeros(size(ExpData.metabolites{3}.conc));
% % % %     end
% % % % 
% % % % 
% % % %     % plotting
% % % %     subplot(7,6,i)
% % % % %         h1.Color = 'black';
% % % %     %plot(tempFF01time,tempFF01conc,'.')
% % % % 
% % % %     nSims = length(selResCell);
% % % %     if nSims == 1
% % % %         colArray = [0 0 0];
% % % %     else
% % % %         colArray = cool(nSims);
% % % %     end
% % % %     for j = 1:nSims
% % % %         T_FF01 = selResCell{j}.T_FF01;
% % % %         Y_FF01 = selResCell{j}.Y_FF01;
% % % %         if j == 1
% % % %             plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'LineWidth', 2)
% % % %         else
% % % %             plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:))
% % % %         end
% % % %         hold on
% % % %     end
% % % %     h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
% % % % 
% % % % 
% % % % %         if i == 6
% % % % %             plot(T_FF01, Y_FF01(:,i), 'Color', 'black')   
% % % % %         else
% % % % %             plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'DisplayName', 'sim') 
% % % % %         end
% % % % 
% % % % %         [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
% % % % %         tempSave1{i}.T_FF01 = T_FF01;
% % % % %         tempSave1{i}.Y_FF01 = Y_FF01;
% % % % %         tempSave1{i}.V_FF01 = V_FF01;
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %     Black line at end of feed
% % % %     plot([20 20], ylim,'k')
% % % %     hold on
% % % %     if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
% % % %         xlim([0 400])
% % % %     elseif(setup.GPdataset.GP1800WT == 1)
% % % %         xlim([0 1800])
% % % %     end
% % % %     title(metNames{i})
% % % % %         if ~isempty(a1)
% % % % %             xtxt = 0.5*max(xlim);
% % % % %             ytxt = 0.8*max(ylim);
% % % % %             text(xtxt,ytxt, ['r^2 = ', num2str(R2, '%3.2f')])
% % % % %             xtxt = 0.7*max(xlim);
% % % % %             ytxt = 0.4*max(ylim);
% % % % %             text(xtxt,ytxt, ['R = ', num2str(R, '%3.2f')])
% % % % %         end
% % % % 
% % % % %     figco = figure(20+i);
% % % % %     h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std); hold on
% % % % %     plot(T_FF01, Y_FF01(:,i), 'Color', 'black') 
% % % % %     title(metNames{i})
% % % % %     xlabel('time'); ylabel('C(mM)')
% % % % %     titl = ['metConc', num2str(i),'.jpg']; 
% % % % %     
% % % % %     
% % % % %     % % % % saveas(figco, titl)
% % % % %     close(figco)
% % % % end
% % % % suptitle('Metabolite concentrations - last cycle');
% % % % % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});
% % % % 
% % % % 
% % % % %%
% % % % 
% % % % if isfield(setup,'glycSynthDeg_separate')
% % % %     if setup.glycSynthDeg_separate == 1
% % % %         % plotting
% % % %         subplot(7,6,40)
% % % % %         h1.Color = 'black';
% % % %         %plot(tempFF01time,tempFF01conc,'.')
% % % % 
% % % %         nSims = length(selResCell);
% % % %         colArray = cool(nSims);
% % % % 
% % % %         for j = 1:nSims
% % % %             T_FF01 = selResCell{j}.T_FF01;
% % % %             Y_FF01 = selResCell{j}.Y_FF01;
% % % %             plot(T_FF01, Y_FF01(:,40), 'Color', colArray(j,:)) 
% % % %             hold on
% % % %         end
% % % %         h1 = plot([0 399],[100 100], ':.b','MarkerSize',10);
% % % % 
% % % %     %     Black line at end of feed
% % % %         plot([20 20], ylim,'k')
% % % %         hold on
% % % % %             xlim([0 400])
% % % %         if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
% % % %             xlim([0 400])
% % % %         elseif(setup.GPdataset.GP1800WT == 1)
% % % %             xlim([0 1800])
% % % %         end
% % % %         title('glycogen')
% % % %     end
% % % % else
% % % % end
% % % % 
% % % % fig = gcf;
% % % % set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% % % % % % % % saveas(fig,'FFLastCycle_mets.jpg')
% % % % % close
% % % % 
% % % % if isfield(setup,'glycSynthDeg')
% % % %     if setup.glycSynthDeg == 1
% % % %         numFlu = 49;
% % % %         legenda.fluxes{49} = 'glycSynthDeg';
% % % %         ExpData.fluxes{49}.time = ExpData.fluxes{48}.time;
% % % %     else
% % % %         numFlu = 48;
% % % %     end
% % % % else
% % % %     numFlu = 48;
% % % % end
% % % % 
% % % % if isfield(setup,'glycSynthDeg_separate')
% % % %     if setup.glycSynthDeg_separate == 1
% % % %         numFlu = 51;
% % % %         legenda.fluxes{49} = 'glycSynthDeg';
% % % %         legenda.fluxes{50} = 'glycSynth';
% % % %         legenda.fluxes{51} = 'glycDeg';
% % % %         ExpData.fluxes{49}.time = ExpData.fluxes{48}.time;
% % % %         ExpData.fluxes{50}.time = ExpData.fluxes{48}.time;
% % % %         ExpData.fluxes{51}.time = ExpData.fluxes{48}.time;
% % % % 
% % % %     end
% % % % else
% % % % end
% % % % %     legend(legenda.metabolites)
% % % % % % % %     legend(namesHits1)
% % % % 
% % % % % %%
% % % % figure%(2002)
% % % % for i = 1:numFlu
% % % % %     for i = 1:48
% % % % 
% % % %     bl = ExpData.fluxes{i};
% % % % %         t_sim = T_FF01;
% % % % %         v_sim = V_FF01(:,i);
% % % % 
% % % % 
% % % %     if isfield(bl, 'rate')
% % % %         expFF01time = bl.time;
% % % %         expFF01flux = bl.rate;
% % % % %             vip = expFF01flux(expFF01time<400);
% % % % %             tip = expFF01time(expFF01time<400);
% % % % %     %         v_r2_exp = interp1(tip, vip, t_sim, 'pchip','extrap');
% % % % %     %         SSE = sum((v_r2_exp - v_sim).^2);
% % % % %     %         SST = sum((v_r2_exp - mean(v_r2_exp)).^2);
% % % % %             v_r2 = interp1(t_sim, v_sim, tip);
% % % % %             SSE = sum((vip - v_r2).^2);
% % % % %             SST = sum((vip - mean(vip)).^2);
% % % % %             R2 = 1 - SSE/SST;
% % % % %             R_matrix = corrcoef(v_r2, vip);
% % % % %     %         R_matrix = corrcoef(v_sim, v_r2_exp);
% % % % %             R = R_matrix(1,2); 
% % % % 
% % % %     else
% % % %         if((isfield(setup.GPdataset,'GP400M'))&&(setup.GPdataset.GP400M == 1))
% % % %             expFF01time = zeros(5,1);
% % % %             expFF01flux = zeros(5,1);
% % % %         else
% % % %             expFF01time = zeros(size(ExpData.fluxes{3}.time));
% % % %             expFF01flux = zeros(size(ExpData.fluxes{3}.rate));
% % % %         end
% % % %     end
% % % % 
% % % % 
% % % %     subplot(7,8,i)
% % % % %         hold on;
% % % % %         plot(T_FF01, V_FF01(:,i), 'Color', 'black')
% % % %     nSims = length(selResCell);
% % % % 
% % % %     if nSims == 1
% % % %         colArray = [0 0 0];
% % % %     else
% % % %         colArray = cool(nSims);
% % % %     end
% % % %     for j = 1:nSims
% % % %         T_FF01 = selResCell{j}.T_FF01;
% % % %         V_FF01 = selResCell{j}.V_FF01;
% % % % %             plot(T_FF01, V_FF01(:,i), 'Color', colArray(j,:))
% % % %         if j == 1
% % % %             plot(T_FF01, V_FF01(:,i), 'Color', 'black', 'LineWidth', 2)
% % % %         else
% % % %             plot(T_FF01, V_FF01(:,i), 'Color', colArray(j,:))
% % % %         end
% % % %         hold on
% % % %     end  
% % % %     % glycogen experimental data
% % % %     if setup.GPdataset.GP400WT == 1 % 400WT
% % % %         if i == 50
% % % %             expFF01time = dataset.FF01.fluxes_times;
% % % %             expFF01flux = dataset.FF01.fluxes{42}';
% % % %         elseif i == 51
% % % %             expFF01time = dataset.FF01.fluxes_times;
% % % % %             expFF01flux = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}';
% % % %             expFF01flux = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}' + dataset.FF01.fluxes{45}';
% % % %         end
% % % %     elseif setup.GPdataset.GP1800WT == 1 % 1800WT
% % % %         if i == 50
% % % %             expFF01time = dataset.FF03.fluxes_times;
% % % %             expFF01flux = dataset.FF03.fluxes{42}';
% % % %         elseif i == 51
% % % %             expFF01time = dataset.FF03.fluxes_times;
% % % %             expFF01flux = dataset.FF03.fluxes{43}' - dataset.FF03.fluxes{44}' + dataset.FF03.fluxes{45}';
% % % %         end
% % % %     elseif setup.GPdataset.GP400M == 1 % 400M
% % % % %             if i == 50
% % % % %                 expFF01time = dataset.FF04.fluxes_times;
% % % % %                 expFF01flux = dataset.FF04.fluxes{42}';
% % % % %             elseif i == 51
% % % % %                 expFF01time = dataset.FF04.fluxes_times;
% % % % %                 expFF01flux = dataset.FF04.fluxes{43}' - dataset.FF01.fluxes{44}' + dataset.FF01.fluxes{45}';
% % % % %             end
% % % %     end
% % % %     plot(expFF01time, expFF01flux,'b.:','MarkerSize',10)      
% % % % 
% % % % 
% % % %     %     Black line at end of feed
% % % %     plot([20 20], ylim,'k'); hold on;
% % % %     title(legenda.fluxes{i})
% % % % %         xlim([0 400])
% % % %     if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
% % % %         xlim([0 400])
% % % %     elseif(setup.GPdataset.GP1800WT == 1)
% % % %         xlim([0 1800])
% % % %     end
% % % %     if isfield(bl, 'rate')
% % % % %             xtxt = 0.5*max(xlim);
% % % % %             ytxt = 0.8*max(ylim);
% % % % %             text(xtxt,ytxt, ['r^2 = ', num2str(R2, '%3.2f')])
% % % % %             xtxt = 0.7*max(xlim);
% % % % %             ytxt = 0.4*max(ylim);
% % % % %             text(xtxt,ytxt, ['R = ', num2str(R, '%3.2f')])
% % % %     end
% % % % 
% % % % end
% % % % suptitle('Enzymatic fluxes - last cycle');
% % % % % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});
% % % % 
% % % % fig = gcf;
% % % % set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% % % % % % % % saveas(fig,'FFLastCycle_fluxes_exp.jpg')
% % % % % close
% % % % 
% % % % % end



