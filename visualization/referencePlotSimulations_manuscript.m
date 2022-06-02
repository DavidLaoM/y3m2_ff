
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % Simulation Metabolites  % % % % % % % % % % % % 
sel_idxs = 1:38;
    sel_idxs([1,22,23,28,29,30,31,32,33,34,35]) = [];
    sel_idxs2 = [     0     1     2     3     4     5     6,...
                      7     8     9    10    11    12    13,...
                     14    15    16    17    18    19    20,...
                      0     0    21    22    23    24     0,...
                      0     0     0     0     0     0     25,...
                     26    27    28];
    
fh2001 = figure(2001);
% for i = 1:38
for first_plot = 1
    % 
    for i = sel_idxs
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
%             tempFF01time = zeros(size(ExpData.metabolites{3}.time));
%             tempFF01conc = zeros(size(ExpData.metabolites{3}.conc));
%             tempFF01std = zeros(size(ExpData.metabolites{3}.conc));
            tempFF01time = [];
            tempFF01conc = [];
            tempFF01std = [];
        end


        % plotting
%         subplot(7,6,i)
%         subplot(4,7,sel_idxs2(i))
        [temp_num1,~] = size(fh2001.Children);
        i2a = temp_num1 + 1;
        subplot(4,7,i2a)
        
        
    %         h1.Color = 'black';
        %plot(tempFF01time,tempFF01conc,'.')

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
%                 plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'LineWidth', 2)
                plot(T_FF01, Y_FF01(:,i), 'Color', [.7 .7 .7], 'LineWidth', 1)
            else
%                 plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:))
                plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'LineWidth', 1)
            end
            hold on
        end
%         h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
        h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, 'ko',...
            'MarkerSize',4,'MarkerFaceColor','w','CapSize',0);
        % Black line at end of feed
        plot([20 20], ylim,'k')
        hold on
        if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
            xlim([0 400])
        elseif(setup.GPdataset.GP1800WT == 1)
            xlim([0 1800])
        end
%         title(metNames{i})
        temp = metNames{i}(1:end-13);
        temp = erase(temp,',');
        title(temp)
    end
    suptitle('Metabolite concentrations - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});
    % 
% % % % %     if isfield(setup,'glycSynthDeg_separate')
% % % % %         if setup.glycSynthDeg_separate == 1
            % plotting
%             subplot(7,6,40)
            subplot(4,7,28)
    %         h1.Color = 'black';
            %plot(tempFF01time,tempFF01conc,'.')

            nSims = length(selResCell);
            colArray = cool(nSims);

            for j = 1:nSims
                T_FF01 = selResCell{j}.T_FF01;
                Y_FF01 = selResCell{j}.Y_FF01;
%                 plot(T_FF01, Y_FF01(:,40), 'Color', colArray(j,:)) 
                if j == 1                
                    plot(T_FF01, Y_FF01(:,40), 'Color', [.7 .7 .7],... 
                        'Color', [.7 .7 .7], 'LineWidth', 1)
                else
                    plot(T_FF01, Y_FF01(:,40), 'Color', 'k',... 
                        'Color', 'k', 'LineWidth', 1)
                end
                hold on
            end
%             h1 = plot([0 399],[100 100], ':.b','MarkerSize',10);
                h1 = plot(0,100, 'ko','MarkerSize',10,...
                    'MarkerSize',4,'MarkerFaceColor','w');

            %     Black line at end of feed
                plot([20 20], ylim,'k')
                hold on
        %             xlim([0 400])
                if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
                    xlim([0 400])
                elseif(setup.GPdataset.GP1800WT == 1)
                    xlim([0 1800])
                end
                ylim([0 110])
                title('glycogen')
% % % % %         end
% % % % %     else
% % % % %     end

    fig = gcf;
    % set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    set(fig, 'Units', 'Normalized', 'Position', [1 0.2370 0.8000 0.6917]);
    % % % % saveas(fig,'FFLastCycle_mets.jpg')
    % close

% % % % %     if isfield(setup,'glycSynthDeg')
% % % % %         if setup.glycSynthDeg == 1
% % % % %             numFlu = 49;
% % % % %             legenda.fluxes{49} = 'glycSynthDeg';
% % % % %             ExpData.fluxes{49}.time = ExpData.fluxes{48}.time;
% % % % %         else
% % % % %             numFlu = 48;
% % % % %         end
% % % % %     else
% % % % %         numFlu = 48;
% % % % %     end
    numFlu = 49;
    legenda.fluxes{49} = 'glycSynthDeg';
    ExpData.fluxes{49}.time = ExpData.fluxes{48}.time;

% % % % %     if isfield(setup,'glycSynthDeg_separate')
% % % % %         if setup.glycSynthDeg_separate == 1
            numFlu = 51;
            legenda.fluxes{49} = 'glycSynthDeg';
            legenda.fluxes{50} = 'glycSynth';
            legenda.fluxes{51} = 'glycDeg';
            ExpData.fluxes{49}.time = ExpData.fluxes{48}.time;
            ExpData.fluxes{50}.time = ExpData.fluxes{48}.time;
            ExpData.fluxes{51}.time = ExpData.fluxes{48}.time;
% % % % %         end
% % % % %     else
% % % % %     end

end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % Simulation Rates  % % % % % % % % % % % % % % % 
sel_idxs = 1:numFlu;
    sel_idxs([16,22,25,29,30,31,32,33,43,44,49]) = [];
fh2002 = figure(2002);
% for i = 1:numFlu
for i = sel_idxs
%     for i = 1:48
    bl = ExpData.fluxes{i};

    % 
    if isfield(bl, 'rate')
        expFF01time = bl.time;
        expFF01flux = bl.rate;
    else
        if((isfield(setup.GPdataset,'GP400M'))&&(setup.GPdataset.GP400M == 1))
            expFF01time = zeros(5,1);
            expFF01flux = zeros(5,1);
        else
            expFF01time = [];
            expFF01flux = [];
        end
    end


%     subplot(7,8,i)
    [temp_num,~] = size(fh2002.Children);
    i2 = temp_num + 1;
    subplot(5,8,i2)
    
    
    nSims = length(selResCell);

    if nSims == 1
        colArray = [0 0 0];
    else
        colArray = cool(nSims);
    end
    for j = 1:nSims
        T_FF01 = selResCell{j}.T_FF01;
        V_FF01 = selResCell{j}.V_FF01;
%             plot(T_FF01, V_FF01(:,i), 'Color', colArray(j,:))
        if j == 1
            plot(T_FF01, V_FF01(:,i), 'Color', [.7 .7 .7], 'LineWidth', 1)
        else
            plot(T_FF01, V_FF01(:,i), 'Color', 'black', 'LineWidth', 1)
        end
        hold on
    end  
    % glycogen experimental data
    if setup.GPdataset.GP400WT == 1 % 400WT
        if i == 50
            expFF01time = dataset.FF01.fluxes_times;
            expFF01flux = dataset.FF01.fluxes{42}';
        elseif i == 51
            expFF01time = dataset.FF01.fluxes_times;
%             expFF01flux = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}';
            expFF01flux = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}' + dataset.FF01.fluxes{45}';
        end
    elseif setup.GPdataset.GP1800WT == 1 % 1800WT
        if i == 50
            expFF01time = dataset.FF03.fluxes_times;
            expFF01flux = dataset.FF03.fluxes{42}';
        elseif i == 51
            expFF01time = dataset.FF03.fluxes_times;
            expFF01flux = dataset.FF03.fluxes{43}' - dataset.FF03.fluxes{44}' + dataset.FF03.fluxes{45}';
        end
    elseif setup.GPdataset.GP400M == 1 % 400M
%             if i == 50
%                 expFF01time = dataset.FF04.fluxes_times;
%                 expFF01flux = dataset.FF04.fluxes{42}';
%             elseif i == 51
%                 expFF01time = dataset.FF04.fluxes_times;
%                 expFF01flux = dataset.FF04.fluxes{43}' - dataset.FF01.fluxes{44}' + dataset.FF01.fluxes{45}';
%             end
    end
    plot(expFF01time, expFF01flux,'ko','MarkerSize',4,...
        'MarkerFaceColor','w','LineWidth',1)      


    %     Black line at end of feed
    plot([20 20], ylim,'k'); hold on;
    % 
%     title(legenda.fluxes{i})
%         xlim([0 400])
    if i2 <= 39
        temp = legenda.fluxes{i}(1:end-20);
        temp = erase(temp,',');
        title(temp)
    else
        title(legenda.fluxes{i})
    end

    if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
        xlim([0 400])
    elseif(setup.GPdataset.GP1800WT == 1)
        xlim([0 1800])
    end
    if isfield(bl, 'rate')
%             xtxt = 0.5*max(xlim);
%             ytxt = 0.8*max(ylim);
%             text(xtxt,ytxt, ['r^2 = ', num2str(R2, '%3.2f')])
%             xtxt = 0.7*max(xlim);
%             ytxt = 0.4*max(ylim);
%             text(xtxt,ytxt, ['R = ', num2str(R, '%3.2f')])
    end

end
suptitle('Enzymatic fluxes - last cycle');
% suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});

fig = gcf;
% set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
set(fig, 'Units', 'Normalized', 'Position', [1 0.2370 0.8000 0.6917]);
% % % % saveas(fig,'FFLastCycle_fluxes_exp.jpg')


%% later edists y-limits
fh2001.Children(24).YLim = [1.0500 1.1500];
fh2001.Children(20).YLim = [0 5];
fh2002.Children(26).YLim = [0.0000 0.0500];
fh2002.Children(21).YLim = [0.0000 0.5000];
fh2002.Children(13).YLim = [0.0007876 0.000788];
fh2002.Children(4).YLim = [-1e-0 0]; %[-6.417e-05 -6.415e-05];

