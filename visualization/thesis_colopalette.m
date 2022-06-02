% % color blind inclusive

% zesty color palette
cpal_small.orange = [245,121,58]/255;%#F5793A % orange
cpal_small.magenta = [169,90,161]/255;%#A95AA1 % magenta
cpal_small.lightblue = [133,192,249]/255;%#85C0F9 % lightblue
cpal_small.darkblue = [15,32,128]/255;%#0F2080 % darkblue

% extended color palette
cpal_broad.black = [0 0 0]/255; % #9F4D3C % black
cpal_broad.orange = [252, 147, 3]/255; % #fc9309 % orange
cpal_broad.skyblue = [38, 165, 203]/255; %#87CEEB % skyblue
% cpal_broad.bluishgreen = [17,102,68]/255; %#116644 % bluishgreen
cpal_broad.bluishgreen = [21,128,85]/255; %#158055 % bluishgreen
% cpal_broad.yellow = [249 198 5]/255; %#f9c605 % yellow
cpal_broad.yellow = [252 233 3]/255; %#fce903 % yellow
% cpal_broad.blue = [0 133 202]/255; %#0085CA % blue
cpal_broad.blue = [0 128 255]/255; %#0080ff % blue
% cpal_broad.vermillion = [227, 66, 52]/255; % #9F4D3C % vermillion
cpal_broad.vermillion = [255, 51, 52]/255; % #ff3300 % vermillion
cpal_broad.reddishpurple = [255 153 204]/255; %#FF94C6 % reddishpurple

%%
temp_visualize = 0;
if temp_visualize == 1
    figure(11),

    % for i = 1:12
    %     subplot(3,4,i)
    %     area([0 1],[1 1],'FaceColor','red')
    % end

    % 
    fn = fieldnames(cpal_small);
    C_small = struct2cell(cpal_small);
    for k=1:numel(fn)
        if( isnumeric(cpal_small.(fn{k})) )
            subplot(3,4,k)
            area([0 1],[1 1],'FaceColor',C_small{k})
        end
    end

    % 
    fn = fieldnames(cpal_broad);
    C_broad = struct2cell(cpal_broad);
    for k=1:numel(fn)
        if( isnumeric(cpal_broad.(fn{k})) )
            subplot(3,4,k+4)
            area([0 1],[1 1],'FaceColor',C_broad{k})
        end
    end
    
end
clear temp_visualize
