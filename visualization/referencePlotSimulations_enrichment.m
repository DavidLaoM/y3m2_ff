% % setup options
% plotMode = 1; % single simulation
% plotMode = 2; % multiple simulations on top of each other (scatter or psa-like)
% plotMode = 10; % all

% Simulations with experimental data % % % % % % % % % % % % % % % % % % % 
if((plotMode == 1)||(plotMode == 10)||(plotMode == 11))
        % %%
    figure(101)
    for i = 1:38

        %FF01

        a1 = ExpData.metabolites{i};
        t_sim = T_FF01;
        c_sim = Y_FF01(:,i);
        if ~isempty(a1)
            tempFF01time = a1.time;
            tempFF01conc = a1.conc;
            cip = tempFF01conc(tempFF01time<400);
            tip = tempFF01time(tempFF01time<400);
            if isfield(a1, 'stdev')
                tempFF01std = a1.stdev;
            else 
                tempFF01std = zeros(size(tempFF01time)); 
            end
            c_r2 = interp1(t_sim, c_sim, tip);
            SSE = sum((cip - c_r2).^2);
            SST = sum((cip - mean(cip)).^2);
            R2 = 1 - SSE/SST;
            R_matrix = corrcoef(c_r2, cip);
            R = R_matrix(1,2);
        else
            tempFF01time = zeros(size(ExpData.metabolites{3}.time));
            tempFF01conc = zeros(size(ExpData.metabolites{3}.conc));
            tempFF01std = zeros(size(ExpData.metabolites{3}.conc));
        end


        % plotting
        subplot(7,6,i)
        h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
%         h1.Color = 'black';
        %plot(tempFF01time,tempFF01conc,'.')
        hold on
        if i == 6
            plot(T_FF01, Y_FF01(:,i), 'Color', 'black')   
        else
            plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'DisplayName', 'sim') 
        end
    %     Black line at end of feed
        plot([20 20], ylim,'k')
        hold on
        xlim([0 400])
        title(metNames{i})
%         if ~isempty(a1)
%             xtxt = 0.5*max(xlim);
%             ytxt = 0.8*max(ylim);
%             text(xtxt,ytxt, ['r^2 = ', num2str(R2, '%3.2f')])
%             xtxt = 0.7*max(xlim);
%             ytxt = 0.4*max(ylim);
%             text(xtxt,ytxt, ['R = ', num2str(R, '%3.2f')])
%         end

    %     figco = figure(20+i);
    %     h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std); hold on
    %     plot(T_FF01, Y_FF01(:,i), 'Color', 'black') 
    %     title(metNames{i})
    %     xlabel('time'); ylabel('C(mM)')
    %     titl = ['metConc', num2str(i),'.jpg']; 
    %     
    %     
    %     % % % % saveas(figco, titl)
    %     close(figco)
    end
    suptitle('Metabolite concentrations - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});

    fig = gcf;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % % % % saveas(fig,'FFLastCycle_mets.jpg')
    % close

% % % % %     if isfield(setup,'glycSynthDeg_separate')
% % % % %         if setup.glycSynthDeg_separate == 1
            % plotting
            subplot(7,6,40)
    %         h1.Color = 'black';
            %plot(tempFF01time,tempFF01conc,'.')

            nSims = length(selResCell);
            colArray = cool(nSims);

            for j = 1:nSims
                T_FF01_b = selResCell{j}.T_FF01;
                Y_FF01_b = selResCell{j}.Y_FF01;
                plot(T_FF01_b, Y_FF01_b(:,40), 'Color', colArray(j,:)) 
            hold on
            end
            h1 = plot([0 399],[100 100], ':.b','MarkerSize',10);

        %     Black line at end of feed
            plot([20 20], ylim,'k')
            hold on
            xlim([0 400])
            title('glycogen')
% % % % %         end
% % % % %     else
% % % % %     end
    
    
    % %%
    figure(102)
    for i = 1:48

        bl = ExpData.fluxes{i};
        t_sim = T_FF01;
        v_sim = V_FF01(:,i);


        if isfield(bl, 'rate')
            expFF01time = bl.time;
            expFF01flux = bl.rate;
            vip = expFF01flux(expFF01time<400);
            tip = expFF01time(expFF01time<400);
    %         v_r2_exp = interp1(tip, vip, t_sim, 'pchip','extrap');
    %         SSE = sum((v_r2_exp - v_sim).^2);
    %         SST = sum((v_r2_exp - mean(v_r2_exp)).^2);
            v_r2 = interp1(t_sim, v_sim, tip);
            SSE = sum((vip - v_r2).^2);
            SST = sum((vip - mean(vip)).^2);
            R2 = 1 - SSE/SST;
            R_matrix = corrcoef(v_r2, vip);
    %         R_matrix = corrcoef(v_sim, v_r2_exp);
            R = R_matrix(1,2); 

        else
            expFF01time = zeros(size(ExpData.fluxes{3}.time));
            expFF01flux = zeros(size(ExpData.fluxes{3}.rate));
        end


        subplot(7,7,i)
        plot(expFF01time, expFF01flux,'b.:','MarkerSize',10)
        hold on;
        plot(T_FF01, V_FF01(:,i), 'Color', 'black')
        %     Black line at end of feed
        plot([20 20], ylim,'k'); hold on;
        title(legenda.fluxes{i})
        xlim([0 400])
        if isfield(bl, 'rate')
            xtxt = 0.5*max(xlim);
            ytxt = 0.8*max(ylim);
            text(xtxt,ytxt, ['r^2 = ', num2str(R2, '%3.2f')])
            xtxt = 0.7*max(xlim);
            ytxt = 0.4*max(ylim);
            text(xtxt,ytxt, ['R = ', num2str(R, '%3.2f')])
        end

    end
    suptitle('Enzymatic fluxes - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});

    fig = gcf;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % % % % saveas(fig,'FFLastCycle_fluxes_exp.jpg')
    % close
end

% Only simulations % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if(plotMode == 10)

    figure(1001)
    for i = 1:38
        subplot(7,6,i)
        plot(T_FF01, Y_FF01(:,i), 'Color', 'black')
        title(metNames{i})
        xlim([0 400])
    end
    suptitle('Metabolite concentrations - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});

    fig = gcf;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % % % % saveas(fig,'FFLastCycle_fluxes.jpg')
    % close

    figure(1002)
    for i = 1:48
        subplot(7,7,i)
        plot(T_FF01, V_FF01(:,i), 'Color', 'black')
        title(legenda.fluxes{i})
        xlim([0 400])
    end
    suptitle('Enzymatic fluxes - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});

    fig = gcf;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % % % % saveas(fig,'FFLastCycle_fluxes.jpg')
    % close

end



% Simulations with exprimental?tracing? data % % % % % % % % % % % % % % % % % % % 
if(plotMode == 2)
        % %%
    figure%(2001)
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
        subplot(7,6,i)
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
                plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'LineWidth', 2)
            else
                plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:))
            end
            hold on
        end
        h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
        
        
%         if i == 6
%             plot(T_FF01, Y_FF01(:,i), 'Color', 'black')   
%         else
%             plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'DisplayName', 'sim') 
%         end
        
%         [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%         tempSave1{i}.T_FF01 = T_FF01;
%         tempSave1{i}.Y_FF01 = Y_FF01;
%         tempSave1{i}.V_FF01 = V_FF01;
        
        
        
        
    %     Black line at end of feed
        plot([20 20], ylim,'k')
        hold on
        if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
            xlim([0 400])
        elseif(setup.GPdataset.GP1800WT == 1)
            xlim([0 1800])
        end
        title(metNames{i})
%         if ~isempty(a1)
%             xtxt = 0.5*max(xlim);
%             ytxt = 0.8*max(ylim);
%             text(xtxt,ytxt, ['r^2 = ', num2str(R2, '%3.2f')])
%             xtxt = 0.7*max(xlim);
%             ytxt = 0.4*max(ylim);
%             text(xtxt,ytxt, ['R = ', num2str(R, '%3.2f')])
%         end

    %     figco = figure(20+i);
    %     h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std); hold on
    %     plot(T_FF01, Y_FF01(:,i), 'Color', 'black') 
    %     title(metNames{i})
    %     xlabel('time'); ylabel('C(mM)')
    %     titl = ['metConc', num2str(i),'.jpg']; 
    %     
    %     
    %     % % % % saveas(figco, titl)
    %     close(figco)
    end
    suptitle('Metabolite concentrations - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});


% % % % %     if isfield(setup,'glycSynthDeg_separate')
% % % % %         if setup.glycSynthDeg_separate == 1
            % plotting
            subplot(7,6,40)
    %         h1.Color = 'black';
            %plot(tempFF01time,tempFF01conc,'.')

            nSims = length(selResCell);
            colArray = cool(nSims);

            for j = 1:nSims
                T_FF01 = selResCell{j}.T_FF01;
                Y_FF01 = selResCell{j}.Y_FF01;
                plot(T_FF01, Y_FF01(:,40), 'Color', colArray(j,:)) 
                hold on
            end
            h1 = plot([0 399],[100 100], ':.b','MarkerSize',10);

        %     Black line at end of feed
            plot([20 20], ylim,'k')
            hold on
%             xlim([0 400])
            if((setup.GPdataset.GP400WT == 1)||(setup.GPdataset.GP400M == 1))
                xlim([0 400])
            elseif(setup.GPdataset.GP1800WT == 1)
                xlim([0 1800])
            end
            title('glycogen')
% % % % %         end
% % % % %     else
% % % % %     end
    
    fig = gcf;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
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
%     legend(legenda.metabolites)
% % % %     legend(namesHits1)
    
    % %%
    figure%(2002)
    for i = 1:numFlu
%     for i = 1:48

        bl = ExpData.fluxes{i};
%         t_sim = T_FF01;
%         v_sim = V_FF01(:,i);


        if isfield(bl, 'rate')
            expFF01time = bl.time;
            expFF01flux = bl.rate;
%             vip = expFF01flux(expFF01time<400);
%             tip = expFF01time(expFF01time<400);
%     %         v_r2_exp = interp1(tip, vip, t_sim, 'pchip','extrap');
%     %         SSE = sum((v_r2_exp - v_sim).^2);
%     %         SST = sum((v_r2_exp - mean(v_r2_exp)).^2);
%             v_r2 = interp1(t_sim, v_sim, tip);
%             SSE = sum((vip - v_r2).^2);
%             SST = sum((vip - mean(vip)).^2);
%             R2 = 1 - SSE/SST;
%             R_matrix = corrcoef(v_r2, vip);
%     %         R_matrix = corrcoef(v_sim, v_r2_exp);
%             R = R_matrix(1,2); 

        else
            if((isfield(setup.GPdataset,'GP400M'))&&(setup.GPdataset.GP400M == 1))
                expFF01time = zeros(5,1);
                expFF01flux = zeros(5,1);
            else
                expFF01time = zeros(size(ExpData.fluxes{3}.time));
                expFF01flux = zeros(size(ExpData.fluxes{3}.rate));
            end
        end


        subplot(7,8,i)
%         hold on;
%         plot(T_FF01, V_FF01(:,i), 'Color', 'black')
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
                plot(T_FF01, V_FF01(:,i), 'Color', 'black', 'LineWidth', 2)
            else
                plot(T_FF01, V_FF01(:,i), 'Color', colArray(j,:))
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
        plot(expFF01time, expFF01flux,'b.:','MarkerSize',10)      
        
        
        %     Black line at end of feed
        plot([20 20], ylim,'k'); hold on;
        title(legenda.fluxes{i})
%         xlim([0 400])
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
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % % % % saveas(fig,'FFLastCycle_fluxes_exp.jpg')
    % close
%     legend(legenda.fluxes)
% % % %     legend(namesHits1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if((isfield(setup.GPdataset,'GP400WT'))&&(setup.GPdataset.GP400WT == 1))
    % %%
    figure%(2003)
    for i = 41:67

        %FF01

        a1 = ExpData.metabolites{i};
%         t_sim = T_FF01;
%         c_sim = Y_FF01(:,i);
        if(~isempty(a1)&&(i ~= 67))
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
        subplot(7,6,i-40)
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
%             plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:)) 
            if j == 1
                plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'LineWidth', 2)
            else
                plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:)) 
            end
            hold on
        end
        if i == 65
            h1 = errorbar(a1.conc_time,tempFF01conc,tempFF01std(1:end-2), ':.r','MarkerSize',10, 'DisplayName', 'exp');
        else
            h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.r','MarkerSize',10, 'DisplayName', 'exp');
        end
        
        
%         if i == 6
%             plot(T_FF01, Y_FF01(:,i), 'Color', 'black')   
%         else
%             plot(T_FF01, Y_FF01(:,i), 'Color', 'black', 'DisplayName', 'sim') 
%         end
        
%         [T_FF01,Y_FF01,V_FF01] = simulate_FF(xSel,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
%         tempSave1{i}.T_FF01 = T_FF01;
%         tempSave1{i}.Y_FF01 = Y_FF01;
%         tempSave1{i}.V_FF01 = V_FF01;
        
        
        
        
    %     Black line at end of feed
        plot([20 20], ylim,'k')
        hold on
        xlim([0 400])
        title(metNames{i})
%         if ~isempty(a1)
%             xtxt = 0.5*max(xlim);
%             ytxt = 0.8*max(ylim);
%             text(xtxt,ytxt, ['r^2 = ', num2str(R2, '%3.2f')])
%             xtxt = 0.7*max(xlim);
%             ytxt = 0.4*max(ylim);
%             text(xtxt,ytxt, ['R = ', num2str(R, '%3.2f')])
%         end

    %     figco = figure(20+i);
    %     h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std); hold on
    %     plot(T_FF01, Y_FF01(:,i), 'Color', 'black') 
    %     title(metNames{i})
    %     xlabel('time'); ylabel('C(mM)')
    %     titl = ['metConc', num2str(i),'.jpg']; 
    %     
    %     
    %     % % % % saveas(figco, titl)
    %     close(figco)
    end
    suptitle('13C LABELLED Metabolite concentrations - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});


%     if isfield(setup,'glycSynthDeg_separate')
%         if setup.glycSynthDeg_separate == 1
%             % plotting
%             subplot(7,6,40)
%     %         h1.Color = 'black';
%             %plot(tempFF01time,tempFF01conc,'.')
% 
%             nSims = length(selResCell);
%             colArray = cool(nSims);
% 
%             for j = 1:nSims
%                 T_FF01 = selResCell{j}.T_FF01;
%                 Y_FF01 = selResCell{j}.Y_FF01;
%                 plot(T_FF01, Y_FF01(:,40), 'Color', colArray(j,:)) 
%             hold on
%             end
%             h1 = plot([0 399],[100 100], ':.b','MarkerSize',10);
% 
%         %     Black line at end of feed
%             plot([20 20], ylim,'k')
%             hold on
%             xlim([0 400])
%             title('glycogen')
%         end
%     else
%     end
    
    fig = gcf;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
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
%     legend(legenda.metabolites)
% % % %     legend(namesHits1)
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % %     % labeling data
% % % %     figure(2003)
% % % %     for i = 41:67
% % % % 
% % % %         %FF01
% % % % %         if i == 67
% % % % %             disp('time to stop')
% % % % %         end
% % % %         a1 = ExpData.metabolites{i};
% % % %         if i == 67
% % % %             ExpData.metabolites{i} = [];
% % % %             a1 = ExpData.metabolites{i};
% % % %         end
% % % % %         t_sim = T_FF01;
% % % % %         c_sim = Y_FF01(:,i);
% % % %         if ~isempty(a1)
% % % %             tempFF01time = a1.time;
% % % %             tempFF01conc = a1.conc;
% % % %             if i == 65
% % % %                 tempFF01time = a1.conc_time;
% % % %             end
% % % % %             cip = tempFF01conc(tempFF01time<400);
% % % % %             tip = tempFF01time(tempFF01time<400);
% % % %             if isfield(a1, 'stdev')
% % % %                 tempFF01std = a1.stdev;
% % % %             else 
% % % %                 tempFF01std = zeros(size(tempFF01time)); 
% % % %             end
% % % % %             c_r2 = interp1(t_sim, c_sim, tip);
% % % % %             SSE = sum((cip - c_r2).^2);
% % % % %             SST = sum((cip - mean(cip)).^2);
% % % % %             R2 = 1 - SSE/SST;
% % % % %             R_matrix = corrcoef(c_r2, cip);
% % % % %             R = R_matrix(1,2);
% % % %         else
% % % %             tempFF01time = zeros(size(ExpData.metabolites{3}.time));
% % % %             tempFF01conc = zeros(size(ExpData.metabolites{3}.conc));
% % % %             tempFF01std = zeros(size(ExpData.metabolites{3}.conc));
% % % %         end
% % % % 
% % % %         % plotting
% % % % %         subplot(7,6,i)
% % % %         i2 = i-40;
% % % %         subplot(7,6,i2)
% % % % %         h1.Color = 'black';
% % % %         %plot(tempFF01time,tempFF01conc,'.')
% % % %         
% % % %         nSims = length(selResCell);
% % % %         if nSims == 1
% % % %             colArray = [0 0 0];
% % % %         else
% % % %             colArray = cool(nSims);
% % % %         end
% % % %         for j = 1:nSims
% % % %             T_FF01 = selResCell{j}.T_FF01;
% % % %             Y_FF01 = selResCell{j}.Y_FF01;
% % % %             plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:)) 
% % % %         hold on
% % % %         end
% % % % %         if i >= 61
% % % % %             disp(i)
% % % % %         end
% % % %         h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
% % % % 
% % % %     %     Black line at end of feed
% % % %         plot([20 20], ylim,'k')
% % % %         hold on
% % % %         xlim([0 400])
% % % %         title(metNames{i})
% % % %     end
% % % %     suptitle('Enrichment metabolite concentrations - last cycle');
% % % %     % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});
% % % %     fig = gcf;
% % % %     set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % labeling data (percentatage)
%     numsPre = [41    42    43    44    45    46    47    48    49,...
%         50    51    52    53    54    55    56    57    58,...
%         59    60    61    62    63    64    65    66    67];
    correspMat = [1    2    3    4    5    6    10    11    12,...
        13    14    17    18    19    20    21    24    25,... %38,...
        26    31    32    33    34    35    36    37    38]; %25];
    correspMat = [1:40, correspMat];

    figure%(2004)
%     figure%(2005)
    for i = 41:67
        %FF01
%         if i == 67
%             disp('time to stop')
%         end
        a1 = ExpData.metabolites{i};
        if i == 67
            ExpData.metabolites{i} = [];
            a1 = ExpData.metabolites{i};
        end
%         t_sim = T_FF01;
%         c_sim = Y_FF01(:,i);
        if ~isempty(a1)
            tempFF01time = a1.time;
            tempFF01conc = a1.conc;
            tempFF01fraction = a1.fraction;
            if i == 65
                tempFF01time = a1.conc_time;
            end
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
%         subplot(7,6,i)
        i2 = i-40;
        subplot(7,6,i2)
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
%             plot(T_FF01, Y_FF01(:,i), 'Color', colArray(j,:)) 
%             plot(T_FF01, Y_FF01(:,i)./Y_FF01(:,correspMat(i)), 'Color', colArray(j,:)) 
%             if i == 59
%                 disp('stop here')
%             end
%             if i == 58
%                 plot(T_FF01, Y_FF01(:,i)./Y_FF01(:,correspMat(i)), 'r-', 'LineWidth', 1.25) 
%             else
%             plot(T_FF01, Y_FF01(:,i)./Y_FF01(:,correspMat(i)), 'r-', 'LineWidth', 1.25)
            
            if j == 1
                plot(T_FF01, Y_FF01(:,i)./Y_FF01(:,correspMat(i)), 'r-', 'LineWidth', 1.25)
            else
                plot(T_FF01, Y_FF01(:,i)./Y_FF01(:,correspMat(i)), '-', 'Color', colArray(j,:))
            end
            
            
%             end
            hold on
        end
%         if i == 58
%             disp(i)
%         end
%         h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
%         disp(i);
%         if sum(tempFF01conc) == 0
%         if length(tempFF01conc) ~= 1
        if((sum(tempFF01conc) == 0)||(length(tempFF01conc) == 1))
%             h1 = errorbar(tempFF01time,tempFF01conc,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
            h1 = plot(tempFF01time, tempFF01conc, 'r+', 'MarkerSize', 5);
        else
%             if i == 66
%                 disp('stop')
%             end
            if i == 65
%                 disp('stop')
                int_val = interp1(ExpData.metabolites{correspMat(i)}.time(1:15), ExpData.metabolites{correspMat(i)}.conc(1:15),tempFF01time,'pchip');
                tempFF01time = datasetEnrich.FF01.metabolites.ECglucose.time;
                tempFF01std = zeros(size(tempFF01time));
            else
                int_val = interp1(ExpData.metabolites{correspMat(i)}.time, ExpData.metabolites{correspMat(i)}.conc,tempFF01time,'pchip');
            end
%             h1 = errorbar(tempFF01time,tempFF01conc./int_val,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
%             h1 = errorbar(tempFF01time,tempFF01fraction,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
            if i ~= 58
                h1 = plot(tempFF01time,tempFF01fraction, 'r+', 'MarkerSize', 5);
            end
        end
        
    %     Black line at end of feed
        plot([20 20], ylim,'k')
        hold on
        xlim([0 400])
        title(metNames{i})
    end
    % addition TREic
    subplot(7,6,28)
    % sims
    temp_t = selResCell{j}.T_FF01;
    temp_y = Y_FF01(:,67) + Y_FF01(:,58);
    temp_y_corresp = Y_FF01(:,38) + Y_FF01(:,25);
    plot(temp_t, temp_y./temp_y_corresp, 'r-', 'LineWidth', 1.25) 
    hold on
    % expData
    plot(dataFF01.metabolites{58}.time,dataFF01.metabolites{58}.fraction, 'r+', 'MarkerSize', 5)
%     errorbar(tempFF01time,tempFF01fraction,tempFF01std, ':.b','MarkerSize',10, 'DisplayName', 'exp');
    title('tre_{ICtotal}')
    
    
    suptitle('(%) Enrichment metabolite concentrations - last cycle');
    % suptitle({'Enzymatic fluxes - last cycle','Opt. cost func. incl. UDPG, G1P, T6P & TRE(cyt+vac) w=1e-2'});
    fig = gcf;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % 
    % 
end

