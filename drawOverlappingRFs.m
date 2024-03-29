function [] = drawOverlappingRFs(h, RF, chArea, xs, ys, tickSp, aList)


% ==== uncomment to make fake gaussian RFs for testing
% chArea = [1 1 1 1 1 1 2 2 2 2 2]; 
% aList= {'V1', 'V2'};
% RF = cell(1, length(chArea)); 
% for ich = 1: length(chArea)
% area = chArea(ich);
% sigRange = [2 5; 5 7];
% xs = -14.5:14.5; ys= -14.5:14.5; 
% [xmat, ymat] = meshgrid(xs, ys); 
% 
% xmu = rand*30-15; ymu = rand*30-15;
% sig = sigRange(area, 1)+rand*diff(sigRange(area,:));
% ypd = makedist('Normal','mu',ymu,'sigma',sig);
% xpd = makedist('Normal','mu',xmu,'sigma',sig);
% 
% RF{ich} = pdf(xpd, xmat).*pdf(ypd, ymat);
% end

figure(h); 
m = length(xs); n = length(ys);

for f = 1:2
    theseChannels = find(chArea == f); %%indices of channels in this area
    nCh = length(theseChannels);

    subplot(1,3,f);    
    
    Array{f} = zeros(n,m); % for each channel
    diffImg = zeros(n,m, 3); % for the composite

    for i = 1:nCh 
        thisRF = RF{theseChannels(i)}.^4;
        if ~isempty(thisRF)
            thisMax = max(thisRF(:));
           binRF = thisRF >= 0.5*thisMax; % threshold the RF 
           if sum(binRF(:)) < numel(binRF)/3 % don't add it if there's too much crap
               Array{f} = Array{f} + double(binRF) ;
           end

        end
    end 

    Array{f} = Array{f}/max(Array{f}(:)); % normalise
    if f==2
        diffImg(:,:,3) = Array{f}; % the third dim of the diffImg is RGB;%need to make it blue and red to be colourblind friendly
        diffImg(:,:,f) = Array{f}; % the third dim of the diffImg is RGB
    else
        diffImg(:,:,f) = Array{f}; % the third dim of the diffImg is RGB
    end
    
    image(xs,ys,(diffImg)); axis image;
     set(gca,'YDir','normal')
    title(aList{f});
    ylabel(['Visual angle (' char(176) ')'])
xlabel(['Visual angle (' char(176) ')'])
%     set(gca, 'XTick', 1:tickSp:length(xs), 'XTickLabel', xs(1:tickSp:end))
%     set(gca, 'YTick', 1:tickSp:length(ys), 'YTickLabel', ys(end:-tickSp:1))

end
    
% finally a subplot that overlays the composites from both areas.
subplot(1,3,3);
diffImg(:,:,1) = Array{1}; % red
diffImg(:,:,2) = Array{2}; % green
diffImg(:,:,3) = Array{2}; % blue - add to make colourblindfriendly
image(xs,ys,(diffImg)); axis image
% set(gca, 'XTick', 1:tickSp:length(xs), 'XTickLabel', xs(1:tickSp:end))
% set(gca, 'YTick', 1:tickSp:length(ys), 'YTickLabel', ys(end:-tickSp:1))
ylabel(['Visual angle (' char(176) ')'])
xlabel(['Visual angle (' char(176) ')'])
 set(gca,'YDir','normal')
title('Overlay');
drawnow 