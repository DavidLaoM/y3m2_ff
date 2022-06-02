% % Control pipeline
%  The aim of this file is to check rerunning codes from start works well
%  after simplification of the pipeline.
check_confirm = 1;
% Please, run this code in sections, as figure numbering is important for
% some plots


%% f1
clearvars -except check_confirm
close all
% tic
F0_reference_simulations
% toc % about 30 secs
if check_confirm == 1
    % 
    fh1_mets_current = figure(2001);
    fh2_rates_current = figure(2002);
    fh4_enrich_current = figure(2004);
    % 
    openfig('F0_FFmets.fig');
%     fh_mets_safecopy = figure(5);
    fh_mets_safecopy = figure(1);
    openfig('F0_FFrates.fig');
%     fh_rates_safecopy = figure(6);
    fh_rates_safecopy = figure(2);
    openfig('F0_FFenrich_mets.fig');
%     fh_enrich_mets_safecopy = figure(7);
    fh_enrich_mets_safecopy = figure(3);
    % 
    if((abs((sum(fh1_mets_current.Children(end).Children(5).YData) - sum(fh_mets_safecopy.Children(end).Children(5).YData))) < 1.25) && (abs((sum(fh2_rates_current.Children(end).Children(4).YData) - sum(fh_rates_safecopy.Children(end).Children(4).YData))) < 1.25) && (abs((sum(fh4_enrich_current.Children(end).Children(4).YData) - sum(fh_enrich_mets_safecopy.Children(end).Children(4).YData))) < 1.25)) 
        disp('F0 reference simulations are reproduced properly.')
    else
        disp('F0 reference simulations are NOT reproduced properly.')
    end
%     close(1:7)
    close(1:3)
end

