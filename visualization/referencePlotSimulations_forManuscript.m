
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % METABOLITE PROFILES % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% clf(2001)
figure(2001)
for i = 1:38

    %FF01

    a1 = ExpData.metabolites{i};
%         t_sim = T_FF01;
%         c_sim = Y_FF01(:,i);
    if ~isempty(a1)
        tempFF01time = a1.time;
        tempFF01conc = a1.conc;
%             cip = tempFF01conc(tempFF01time<400);
%             tip = tempFF01time(tempFF01time<400);
        if isfield(a1, 'stdev')
            tempFF01std = a1.stdev;
        else 
            tempFF01std = zeros(size(tempFF01time)); 
        end
%             c_r2 = interp1(t_sim, c_sim, tip);
%             SSE = sum((cip - c_r2).^2);
%             SST = sum((cip - mean(cip)).^2);
%             R2 = 1 - SSE/SST;
%             R_matrix = corrcoef(c_r2, cip);
%             R = R_matrix(1,2);
    else
        tempFF01time = zeros(size(ExpData.metabolites{3}.time));
        tempFF01conc = zeros(size(ExpData.metabolites{3}.conc));
        tempFF01std = zeros(size(ExpData.metabolites{3}.conc));
    end
    
    % plotting
    subplot(6,7,i)
    % area shaded plot + send to the background
    temp = [max(abs(selResCell{1}.Y_FF01(:,i))), max(abs(tempFF01conc))];
    ymax = max(temp);
    ymax_inc = 1.3;
    area([0 20],[ymax*ymax_inc ymax*ymax_inc], 'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
    hold on   

    % simulated and experimental data
    nSims = length(selResCell);
    if nSims == 1
        colArray = [0 0 0];
    else
        colArray = cool(nSims);
    end
    % 
%     h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, 'o', 'MarkerSize', 3.5, 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'Color','r');
    for j = 1:nSims
        T_FF01 = selResCell{j}.T_FF01;
        Y_FF01 = selResCell{j}.Y_FF01;
        if j == 1
            plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'LineWidth', 1.75)
        else
            plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:))
        end
        hold on
    end
    h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, 'o', 'MarkerSize', 5,...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'Color','k',...
        'LineWidth', 1);
    h1.CapSize = 0.5;
    h1b = plot(tempFF01time,tempFF01conc, 'ko', 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5);
    
    % layout
    xlim([0 400])% set x-limits
    xticks([0 200 400])
    if ymax == 0, ymax = 1; end
    ylim([0 ymax*ymax_inc])% set y-limits
    tempText = metNames{i};
    if i == 14 % GLYCERAL3P to GAP
        tempText = 'GAP';
    elseif i == 24 % UDP_GLC to UDPGLC
        tempText = 'UDPGlc';
    elseif i == 3 % FBP
        tempText = 'FBP';
    elseif i == 6 % GLCic
        tempText = 'GLCic';
    elseif i == 19 % GLYC
        tempText = 'GLYC';
    elseif i == 25 % TREcyt
        tempText = 'TREcyt';
    elseif i == 27 % Pi
        tempText = 'Pi';
    elseif i == 31 % ETOHext
        tempText = 'ETOHext';
    elseif i == 32 % GLYCext
        tempText = 'GLYCext';
    elseif i == 33 % FRUCic
        tempText = 'FRUCic';
    elseif i == 34 % FRUCec
        tempText = 'FRUCec';
    elseif i == 35 % SUCCec
        tempText = 'SUCCec';
    elseif i == 36 % GLCext
        tempText = 'GLCext';
    elseif i == 37 % TREext
        tempText = 'TREext';
    elseif i == 38 % TREvac
        tempText = 'TREvac';
    end
    temp = text(400*0.95, ymax*0.90, strtok(tempText, ','), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    temp2 = text(400*0.95, ymax*0.05, num2str(i), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
end
suptitle('Metabolite concentrations - last cycle');
% suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});
%
% % add GLYCOGEN + TREic + Energy charge
for metaboliteAdditions = 1
    
    % % TREic
    % recall data
    time_sim = selResCell{1}.T_FF01;
    TREic_sim = selResCell{1}.Y_FF01(:,25) + selResCell{1}.Y_FF01(:,38); % TREcyt + TREvac
    time_exp = dataset.FF01.metabolites.ICtreh.time; %ExpData.metabolites{38}.time;
    TREic_exp = dataset.FF01.metabolites.ICtreh.conc; %ExpData.metabolites{25}.conc + ExpData.metabolites{38}.conc; % TREcyt + TREvac
    TREic_stdev = dataset.FF01.metabolites.ICtreh.stdev;
    % plot
    subplot(6,7,39)
    % area
    temp = [max(abs(TREic_exp)), max(abs(TREic_sim))];
    ymax = max(temp);
    area([0 20],[ymax*ymax_inc ymax*ymax_inc], 'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
    % exp and sim data
    hold on
%     plot(time_sim, TREic_sim, 'Color', 'black', 'LineWidth', 2) % simulation
    plot(time_sim, TREic_sim, 'Color', 'black', 'LineWidth', 1.75) % simulation
    % 
    h1 = errorbar(time_exp, TREic_exp, TREic_stdev, 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'k',...
        'LineWidth', 1);
    h1.CapSize = 0.5; % experimental data
    h1b = plot(time_exp, TREic_exp, 'ko', 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5);
%     h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, 'o', 'MarkerSize', 5,...
%         'MarkerFaceColor', 'w', 'MarkerEdgeColor','k', 'Color','k',...
%         'LineWidth', 1);
%     h1.CapSize = 0.5;
%     plot(time_sim, TREic_sim, 'Color', 'black', 'LineWidth', 2) % simulation
    % layout
    xlim([0 400])% set x-limits
    xticks([0 200 400])
    ylim([0 ymax*ymax_inc])% set y-limits
    tempText = 'TREic';
    temp = text(400*0.95, ymax*0.90, strtok(tempText, ','), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    temp2 = text(400*0.95, ymax*0.05, num2str(39), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    
    % % GLYCOGEN
    % plotting
    subplot(6,7,40)
    % 
    area([0 20],[130 130], 'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
    hold on   

    nSims = length(selResCell);
    colArray = cool(nSims);

    for j = 1:nSims
        T_FF01 = selResCell{j}.T_FF01;
        Y_FF01 = selResCell{j}.Y_FF01;
        plot(T_FF01, Y_FF01(:,40), 'Color', 'black', 'LineWidth', 2) 
        hold on
    end
    h1 = errorbar([0 399],[100 100],[0 0], 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'k',...
        'LineWidth', 1);
    h1.CapSize = 0.5;
    h1b = plot([0 399],[100 100], 'ko', 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5);
    
    % layout
    xlim([0 400])% set x-limits
    xticks([0 200 400])
    ylim([0 130])% set y-limits
    temp = text(400*0.95, 100*0.90, 'GLYCO', ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    temp2 = text(400*0.95, ymax*0.05, num2str(40), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    
    % % energy charge = (ATP + ADP/2) / (ATP + ADP + AMP)
    % recall data
    time_sim = selResCell{1}.T_FF01;
    ATP_sim = selResCell{1}.Y_FF01(:,9);
    ADP_sim = selResCell{1}.Y_FF01(:,15);
    AMP_sim = selResCell{1}.Y_FF01(:,16);
    Echarge_sim = (ATP_sim + ADP_sim ./ 2) ./ (ATP_sim + ADP_sim + AMP_sim);
    time_exp = dataset.FF01.metabolites.ICATP.time;
    ATP_exp = dataset.FF01.metabolites.ICATP.conc;
    ADP_exp = dataset.FF01.metabolites.ICADP.conc;
    AMP_exp = dataset.FF01.metabolites.ICAMP.conc;
    Echarge_exp = (ATP_exp + ADP_exp ./ 2) ./ (ATP_exp + ADP_exp + AMP_exp);
    % plot
    subplot(6,7,41)
    % area
    temp = [max(abs(Echarge_exp)), max(abs(Echarge_sim))];
    ymax = max(temp);
    area([0 20],[ymax*ymax_inc ymax*ymax_inc], 'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
    % exp and sim data
    hold on
    plot(time_sim, Echarge_sim, 'Color', 'black', 'LineWidth', 1.75) % simulation
    h1 = errorbar(time_exp,Echarge_exp,zeros(size(Echarge_exp)), 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'k',...
        'LineWidth', 1);
    h1.CapSize = 0.5; % experimental data
    h1b = plot(time_exp, Echarge_exp, 'ko', 'MarkerFaceColor', 'w', ...
        'MarkerSize', 5);
%     plot(time_sim, Echarge_sim, 'Color', 'black', 'LineWidth', 1.75) % simulation
    % layout
    xlim([0 400])% set x-limits
    xticks([0 200 400])
    ylim([0.5 ymax*ymax_inc])% set y-limits
    tempText = 'Energy charge';
    temp = text(400*0.95, 0.15 + ymax*0.90, strtok(tempText, ','), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    temp2 = text(400*0.95, 0.5 + ymax*0.05, num2str(41), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    
end


%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % REACTION RATES PROFILES % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% clf(2002)
figure(2002)
% for i = 1:numFlu
for i = 1:48
    
    % 
    bl = ExpData.fluxes{i};
    if isfield(bl, 'rate')
        expFF01time = bl.time;
        expFF01flux = bl.rate;
    else
        expFF01time = zeros(size(ExpData.fluxes{3}.time));
        expFF01flux = zeros(size(ExpData.fluxes{3}.rate));
    end
    
    % 
    subplot(7,8,i)
    % area shaded plot + send to the background
    temp = [max(abs(selResCell{1}.V_FF01(:,i))), max(abs(expFF01flux))];
    ymax = max(temp);
    ymax_inc = 1.3;
    area([0 20],[ymax*ymax_inc ymax*ymax_inc], 'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
    hold on      
    area([0 20],[-ymax*ymax_inc -ymax*ymax_inc], 'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
   
    % plotting simulated and calculated reaction rates
    nSims = length(selResCell);
    if nSims == 1
        colArray = [0 0 0];
    else
        colArray = cool(nSims);
    end
    % 
    if i == 17 % PGM fwd rev issue
        expFF01time = dataset.FF01.fluxes_times;
        expFF01flux = dataset.FF01.fluxes{39}';
    elseif i == 19 % PGM fwd rev issue
        expFF01time = dataset.FF01.fluxes_times;
        expFF01flux = dataset.FF01.fluxes{12}';
    elseif i == 50
        expFF01time = dataset.FF01.fluxes_times;
        expFF01flux = dataset.FF01.fluxes{42}';
    elseif i == 51
        expFF01time = dataset.FF01.fluxes_times;
        expFF01flux = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}' + dataset.FF01.fluxes{45}';
    end
%     h1 = errorbar(expFF01time,expFF01flux,zeros(size(expFF01flux)), 's', ...
%         'MarkerSize', 6, 'MarkerFaceColor','w', 'MarkerEdgeColor','k', ...
%         'Color','r','LineWidth', 1);
%     h1.CapSize = 0;
    
    for j = 1:nSims
        T_FF01 = selResCell{j}.T_FF01;
        V_FF01 = selResCell{j}.V_FF01;
        if j == 1
            plot(T_FF01, V_FF01(:,i), 'Color', 'black', 'LineWidth', 1.75)
        else
            plot(T_FF01, V_FF01(:,i), 'Color', colArray(j,:))
        end
        hold on
    end  
    % 
    h1 = errorbar(expFF01time,expFF01flux,zeros(size(expFF01flux)), 's', ...
        'MarkerSize', 6, 'MarkerFaceColor','w', 'MarkerEdgeColor','k', ...
        'Color','r','LineWidth', 1.25);
    h1.CapSize = 0;
    
    
    
%     % Black line at end of feed
%     title(legenda.fluxes{i})
%     xlim([0 400])
    
    % layout
    xlim([0 400])% set x-limits
    xticks([0 200 400])
    if ymax == 0, ymax = 1; end
    % 
    if i == 41 % ADH1
        ylim([0 0.46*ymax_inc])% set y-limits
    elseif i == 47 % AGT1
        ylim([-0.001 ymax*ymax_inc])% set y-limits
    elseif i == 48 % 
        ylim([-ymax*ymax_inc 0])% set y-limits
        ymax = -3E-5;
    else
        ylim([0 ymax*ymax_inc])% set y-limits
    end
    % 
    tempText = legenda.fluxes{i};
    tempText = strtok(tempText, ',');
    tempText = strtok(tempText, '}');
    if i == 23
        tempText = 'ETOHtr';
    elseif i == 24
        tempText = 'GLYCEtr';
    elseif i < 49
        tempText = tempText(4:end);
    end
    temp = text(400*0.95, ymax*0.90, tempText, ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
%     % 
%     disp(strtok(tempText, ','))
    %
    if((i == 41)||(i == 47)||(i == 48))
%         e
    else
        temp2 = text(400*0.95, ymax*0.05, num2str(i), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner     
    end
end
suptitle('Enzymatic reaction rates - last cycle');

fig = gcf;
% set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % ENRICHMENT CONCENTRATIONS (%) % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


correspMat = [1    2    3    4    5    6    10    11    12,...
    13    14    17    18    19    20    21    24    25,... %38,...
    26    31    32    33    34    35    36    37    38]; %25];
correspMat = [1:40, correspMat];

% clf(2004)
figure(2004)
%     figure%(2005)
for i = 41:67
    %FF01
    a1 = ExpData.metabolites{i};
    if i == 67
        ExpData.metabolites{i} = [];
        a1 = ExpData.metabolites{i};
    end
    % 
    if ~isempty(a1)
        tempFF01time = a1.time;
        tempFF01conc = a1.conc;
        tempFF01fraction = a1.fraction;
        if i == 65
            tempFF01time = a1.conc_time;
        end
        %
        if isfield(a1, 'stdev')
            tempFF01std = a1.stdev;
        else 
            tempFF01std = zeros(size(tempFF01time)); 
        end
    else
        tempFF01time = zeros(size(ExpData.metabolites{3}.time));
        tempFF01conc = zeros(size(ExpData.metabolites{3}.conc));
        tempFF01std = zeros(size(ExpData.metabolites{3}.conc));
    end

    % plotting
    i2 = i-40;
    subplot(5,6,i2)
    % area shaded plot + send to the background
%     temp = [max(abs(selResCell{1}.V_FF01(:,i))), max(abs(tempFF01conc))];
%     ymax = max(temp);
%     ymax_inc = 1.3;
    area([0 20],[1 1], 'FaceColor', [.9 .9 .9], 'LineStyle', 'none')
    hold on      
    % 
    nSims = length(selResCell);
    if nSims == 1
        colArray = [0 0 0];
    else
        colArray = cool(nSims);
    end
    for j = 1:nSims
        T_FF01 = selResCell{j}.T_FF01;
        Y_FF01 = selResCell{j}.Y_FF01;
        if j == 1
            plot(T_FF01, Y_FF01(:,i)./Y_FF01(:,correspMat(i)), 'k-', 'LineWidth', 1.25)
        else
            plot(T_FF01, Y_FF01(:,i)./Y_FF01(:,correspMat(i)), '-', 'Color', colArray(j,:))
        end
        hold on
    end
    %
    if((sum(tempFF01conc) == 0)||(length(tempFF01conc) == 1))
        h1 = plot(tempFF01time, tempFF01conc, 'k+', 'MarkerSize', 5);
    else
        %
        if i == 65
            int_val = interp1(ExpData.metabolites{correspMat(i)}.time(1:15), ExpData.metabolites{correspMat(i)}.conc(1:15),tempFF01time,'pchip');
            tempFF01time = datasetEnrich.FF01.metabolites.ECglucose.time;
            tempFF01std = zeros(size(tempFF01time));
        else
            int_val = interp1(ExpData.metabolites{correspMat(i)}.time, ExpData.metabolites{correspMat(i)}.conc,tempFF01time,'pchip');
        end
        if i ~= 58
            h1 = plot(tempFF01time,tempFF01fraction, 'k+', 'MarkerSize', 5,...
                'LineWidth', 1.5);
        end
    end

%     Black line at end of feed
%     plot([20 20], ylim,'k')
    hold on
    xlim([0 400])
    ylim([0 1])
%     title(metNames{i})
    if i == 51 % GAP
        tempName = 'GAP';
    elseif i == 53 % GLYC3P
        tempName = 'GLYC3P';
    elseif i == 54 % GLYC
        tempName = 'GLYC';
    elseif i == 57 % UDPGlc
        tempName = 'UDPGlc';
    elseif i == 60 % ETOHec
        tempName = 'ETOHec';
    else
        tempName = strtok(metNames{i}, '#');
    end
%     title(tempName)
    
    temp = text(400*0.95, 0.75, tempName, ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
    temp2 = text(400*0.95, 0.05, num2str(i), ...
        'HorizontalAlignment', 'right', 'FontSize',12, ...
        'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 

end

% addition TREic
subplot(5,6,28)
% sims
temp_t = selResCell{j}.T_FF01;
temp_y = Y_FF01(:,67) + Y_FF01(:,58);
temp_y_corresp = Y_FF01(:,38) + Y_FF01(:,25);
plot(temp_t, temp_y./temp_y_corresp, 'k-', 'LineWidth', 1.25) 
hold on
% expData
plot(dataFF01.metabolites{58}.time,dataFF01.metabolites{58}.fraction, 'k+', 'MarkerSize', 5,...
                'LineWidth', 1.5);
ylim([0 1])
% title('tre_{ICtotal}')
temp = text(400*0.95, 0.75, 'TREic', ...
    'HorizontalAlignment', 'right', 'FontSize',12, ...
    'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 
temp2 = text(400*0.95, 0.05, num2str(68), ...
    'HorizontalAlignment', 'right', 'FontSize',12, ...
    'VerticalAlignment', 'bottom');% label: delete units, x + rather than metabolite name as title, add it as text 0.75 left up corner 

% 
suptitle('(%) Enrichment metabolite concentrations - last cycle');
fig = gcf;
% set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

