function ax = beautifyPlot(fS,ax)
%% This function takes an axis handle and applies standard parameters to it
if nargin<2
    ax = gca;
    fig=gcf;
end
if ~nargin
    fS = 21;
end
%% Standard Settings
fontSize = fS;
lineWidth = 1.9;
figureSize = [.1, .1, .34, .57];%left bottom width height
%% Apply Settings
box off
grid off
set(ax,'FontSize',fontSize,...
       'LineWidth',lineWidth,...
       'TickDir','out');  
    set(fig, 'Units', 'Normalized'); % First change to normalized units.
set(fig, 'OuterPosition', figureSize); % [xLeft, yBottom, width, height]
end