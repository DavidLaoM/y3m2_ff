% Re open simulation figures
clc, clear, close all
openfig('simulations_Y3M2.fig');
% h = get(1,'children');
close(3)
    % figure(1) = concentrations
    % figure(2) = reaction rates
    % figure(4) = enrichment profiles
manus_labels

%% (3000) Enrichment (5) plot
if exist('fh_3000','var')
    clf(3000)
else
    fh_3000 = figure(3000);
end
fh_3000.Position = [1981 284 455 693];

% %% get objects and place them in the right spot
fh3 = get(4,'children');
data2copy3 = ["G6P","F6P","DHAP", ... 
    "P3G","PEP", ...
    "G1P","UDPGlc","T6P","68"];
changePos = cell(1,length(data2copy3));
% 
ax1 = subplot(5,3,1,'parent',fh_3000);
ax2 = subplot(5,3,2,'parent',fh_3000);
ax3 = subplot(5,3,3,'parent',fh_3000);
%     
ax4 = subplot(5,3,4,'parent',fh_3000);
ax5 = subplot(5,3,5,'parent',fh_3000);
ax6 = subplot(5,3,6,'parent',fh_3000);
%     
ax9 = subplot(5,3,9,'parent',fh_3000);
%     
ax12 = subplot(5,3,12,'parent',fh_3000);
%     
ax15 = subplot(5,3,15,'parent',fh_3000);
% 
% loop to search for each metabolite in the name list created
for i = 1:length(data2copy3)
    % loop to screen the entire hf1 object for that metabolite 
    for j = 1:length(fh3)
        sub_fh3 = findobj(fh3(j), 'Type', 'Text');
        % loop to screen in each text field
        for k = 1:length(sub_fh3)
            % ++ if found, copy paste it in figure 3000
            if sub_fh3(k).String == data2copy3(i) % (sub_fh3(k).String == "68")
                    newh = copyobj(fh3(j),3000);
                    % delete x tick axes
% % % %                     newh.XTick = [];
% % % %                     newh.YTick = [];
                    % delete number
                    current_label = findobj(newh, 'Type', 'Text');
                    delete(current_label(1));
                    % give it an adjusted location
                    if i == 1
                        set(newh,'Position',get(ax3,'position'))
                        delete(ax3);
                    elseif i == 2
                        set(newh,'Position',get(ax6,'position'))
                        delete(ax6);
                    elseif i == 3
                        set(newh,'Position',get(ax9,'position'))
                        delete(ax9);
                    elseif i == 4
                        set(newh,'Position',get(ax12,'position'))
                        delete(ax12);
                    elseif i == 5
                        set(newh,'Position',get(ax15,'position'))
                        delete(ax15);
                    elseif i == 6
                        set(newh,'Position',get(ax2,'position'))
                        delete(ax2);
                    elseif i == 7
                        set(newh,'Position',get(ax1,'position'))
                        delete(ax1);
                    elseif i == 8
                        set(newh,'Position',get(ax5,'position'))
                        delete(ax5);
                    elseif i == 9
                        set(newh,'Position',get(ax4,'position'))
                        delete(ax4);
                    end
                    
                    
                    
                    % change line width
                    tempLine = findobj(newh, 'Type', 'Line');
                    tempLine(1).LineWidth = 1;
                    tempLine(2).LineWidth = 2;
            end
        end
    end
end
% %%

%% Saving outcome
% % % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F7A');

