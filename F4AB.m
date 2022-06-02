% Re open simulation figures
% close all
clc, clear, close all
openfig('simulations_Y3M2.fig');
% h = get(1,'children');
close(3)
    % figure(1) = concentrations
    % figure(2) = reaction rates
    % figure(4) = enrichment profiles
manus_labels


%%
% get objects and place them in the right spot
fh1 = get(1,'children');
data2copy = ["G6P","F6P","FBP", ... 
    "GAP","DHAP","GLYC3P", ...
    "BPG","P3G","P2G","PEP", ...
    "Energy charge", ...
    "G1P","UDPGlc","T6P","TREic", ...
    "TREext","GLCext"];

% ax84
if exist('fh_1000','var')
    clf(1000)
end
fh_1000 = figure(1000);
fh_1000.Position = [1371 42 549 954];
% 
ax1 = subplot(8,4,1,'parent',fh_1000);
ax3 = subplot(8,4,3,'parent',fh_1000);
%     
ax5 = subplot(8,4,5,'parent',fh_1000);
ax6 = subplot(8,4,6,'parent',fh_1000);
ax7 = subplot(8,4,7,'parent',fh_1000);
ax8 = subplot(8,4,8,'parent',fh_1000);
% 
ax9 = subplot(8,4,9,'parent',fh_1000);
ax10 = subplot(8,4,10,'parent',fh_1000);
ax11 = subplot(8,4,11,'parent',fh_1000);
% 
ax13 = subplot(8,4,13,'parent',fh_1000);
ax14 = subplot(8,4,14,'parent',fh_1000);
ax15 = subplot(8,4,15,'parent',fh_1000);
ax16 = subplot(8,4,16,'parent',fh_1000);
% 
ax20 = subplot(8,4,20,'parent',fh_1000);
% 
ax24 = subplot(8,4,24,'parent',fh_1000);
% 
ax28 = subplot(8,4,28,'parent',fh_1000);
% 
ax32 = subplot(8,4,32,'parent',fh_1000);



% loop to search for each metabolite in the name list created
for i = 1:length(data2copy) % [1 2 3 4 17 7:10]
    % loop to screen the entire hf1 object for that metabolite 
    for j = 1:length(fh1)
        sub_fh1 = findobj(fh1(j), 'Type', 'Text');
        % loop to screen in each text field
        for k = 1:length(sub_fh1)
            % ++ if found, copy paste it in figure 1000
            if sub_fh1(k).String == data2copy(i)
                newh = copyobj(fh1(j),1000);
                % adjust marker size (one on top of the other)
                newh.Children(4).MarkerSize = 3.5;
                newh.Children(3).MarkerSize = 3.5;
                % delete x tick axes
%                 newh.XTick = [];
                % delete number
                current_label = findobj(newh, 'Type', 'Text');
                delete(current_label(1));
                % change size text
                if i == 6
                    current_label(2).String = 'G3P';
                end
                % set position
                if i == 1
                    set(newh,'Position',get(ax7,'position'))
                    delete(ax7);
                elseif i == 2
                    set(newh,'Position',get(ax11,'position'))
                    delete(ax11);
                elseif i == 3
                    set(newh,'Position',get(ax15,'position'))
                    delete(ax15);
                elseif i == 4
                    set(newh,'Position',get(ax16,'position'))
                    delete(ax16);
                elseif i == 5
                    set(newh,'Position',get(ax14,'position'))
                    delete(ax14);
                elseif i == 6
                    set(newh,'Position',get(ax13,'position'))
                    delete(ax13);
                elseif i == 7
                    set(newh,'Position',get(ax20,'position'))
                    delete(ax20);
                elseif i == 8
                    set(newh,'Position',get(ax24,'position'))
                    delete(ax24);
                elseif i == 9
                    set(newh,'Position',get(ax28,'position'))
                    delete(ax28);
                elseif i == 10
                    set(newh,'Position',get(ax32,'position'))
                    delete(ax32);
                    xlabel(Y3M1_labels.time)
                    ylabel(Y3M1_labels.conc_mMs)
                elseif i == 11
                    set(newh,'Position',get(ax8,'position'))
                    delete(ax8);
                elseif i == 12
                    set(newh,'Position',get(ax10,'position'))
                    delete(ax10);
                elseif i == 13
                    set(newh,'Position',get(ax6,'position'))
                    delete(ax6);
                elseif i == 14
                    set(newh,'Position',get(ax9,'position'))
                    delete(ax9);
                elseif i == 15
                    set(newh,'Position',get(ax5,'position'))
                    delete(ax5);
                elseif i == 16
                    set(newh,'Position',get(ax1,'position'))
                    delete(ax1);
                elseif i == 17
                    set(newh,'Position',get(ax3,'position'))
                    delete(ax3);
                end
                % adjust y-axis 
%                 newh
            end
        end
    end
end

%% saving
% % % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F4ABdata');


%% (1000) Concentrations plot
% open the figure where we want to plot on top (figure 1000)
% % % % I = imread('Y3M2_modelDiagram_simplified_grey.pdf-Page1.png');
% I = imread('Y3M2_modelDiagram_vHeerden2014_alike grey.png');
% clf(1000)
if exist('fh_1001','var')
    clf(1001)
end
fh_1001 = figure(1001);
fh_1001.Position = [1371 42 549 954];

% get objects and place them in the right spot
fh2 = get(2,'children');
data2copy2 = ["GLT","GLK","PGI","PFK","ALD", ... %1 - - %4 -
    "TDH1","PGK","PGM","ENO","PYK","PDC", ... %6 - - - - %11 
    "PGM1","UGP","TPS1","TPS2","NTH1", ... %12 %13 %14 %15 %16
    "glycSynth","glycDeg"]; %17 %18
% 
ax1 = subplot(8,4,1,'parent',fh_1001);
ax2 = subplot(8,4,2,'parent',fh_1001);
ax3 = subplot(8,4,3,'parent',fh_1001);
%     
ax5 = subplot(8,4,5,'parent',fh_1001);
ax6 = subplot(8,4,6,'parent',fh_1001);
ax7 = subplot(8,4,7,'parent',fh_1001);
% 
ax9 = subplot(8,4,9,'parent',fh_1001);

ax11 = subplot(8,4,11,'parent',fh_1001);
% 
ax13 = subplot(8,4,13,'parent',fh_1001);
ax14 = subplot(8,4,14,'parent',fh_1001);
ax15 = subplot(8,4,15,'parent',fh_1001);

% loop to search for each metabolite in the name list created
for i = [1 4, 6 11, 12 13 14 15 16, 17 18] %1:length(data2copy2)
    % loop to screen the entire hf1 object for that metabolite 
    for j = 1:length(fh2)
        sub_fh2 = findobj(fh2(j), 'Type', 'Text');
        % loop to screen in each text field
        for k = 1:length(sub_fh2)
            % ++ if found, copy paste it in figure 1000
            if sub_fh2(k).String == data2copy2(i)
                % 
                newh = copyobj(fh2(j),1001);
                % adjust marker size (one on top of the other)
                newh.Children(4).MarkerSize = 3.5;
                newh.Children(3).MarkerSize = 3.5;
                % delete number
                current_label = findobj(newh, 'Type', 'Text');
                delete(current_label(1));
                % 
                if i == 1
                    set(newh,'Position',get(ax3,'position'))
                    delete(ax3);
                elseif i == 4
                    set(newh,'Position',get(ax7,'position'))
                    delete(ax7);
                elseif i == 6
                    set(newh,'Position',get(ax11,'position'))
                    delete(ax11);
                elseif i == 11
                    set(newh,'Position',get(ax15,'position'))
                    delete(ax15);
                elseif i == 12
                    set(newh,'Position',get(ax6,'position'))
                    delete(ax6);
                elseif i == 13
                    set(newh,'Position',get(ax5,'position'))
                    delete(ax5);
                elseif i == 14
                    set(newh,'Position',get(ax14,'position'))
                    delete(ax14);
                elseif i == 15
                    set(newh,'Position',get(ax13,'position'))
                    delete(ax13);
                elseif i == 16
                    set(newh,'Position',get(ax9,'position'))
                    delete(ax9);
                elseif i == 17
                    set(newh,'Position',get(ax1,'position'))
                    delete(ax1);
                elseif i == 18
                    set(newh,'Position',get(ax2,'position'))
                    delete(ax2);
                end
            end
        end
    end
end

%% saving
% % % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch5_Y3M2_FF\workingVersion\IntraCycles\RESULTS\figures_manuscript\F4Bdata');

