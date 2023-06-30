function [lineOut, fillOut] = stdshade(amatrix,alpha,acolor,F,smth,ax,varargin)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23


if ~exist('acolor','var') || isempty(acolor)
    acolor='r'; 
end

if ~exist('ax','var')
    ax=gca;
end

if ~exist('F','var') || isempty(F)
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1; %no smoothing by default
end  

if ne(size(F,1),1)
    F=F';
end


if smth > 1
    amatrix = boxFilter(amatrix,smth); %use boxfilter to smooth data
end
amean = nanmean(amatrix,1); %get man over first dimension
%astd = nanstd(amatrix,[],1); % to get std shading
astd = nanstd(amatrix,[],1)./sqrt(sum(~isnan(amatrix),1)); % to get sem shading
amean_nan=amean;
std_nan=astd;
astd(isnan(astd))=0;
amean(isnan(amean))=0;
if exist('alpha','var')==0 || isempty(alpha) 
    fillOut = fill(ax,[F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
    acolor='k';
else
    if ~isempty(varargin) && strcmp(varargin{1},'laminar')
        fillOut = fill(ax,[amean+astd fliplr(amean-astd)],[F fliplr(F)],acolor, 'FaceAlpha', alpha,'linestyle','none');%VERTICAL
    else
        fillOut = fill(ax,[F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');%HORIZONTAL
    end
end
fillOut.Annotation.LegendInformation.IconDisplayStyle = 'off';
if ishold==0
    check=true; else check=false;
end

hold on;
% b=ones(3,1)./3;
%     datnonan=amean;
%     endnonan=find(diff(isnan(datnonan))==1);
%     startnonan=find(diff(isnan(datnonan))==-1);
%     [~,p]=max(endnonan-startnonan);
%     lastnonan=endnonan(p);
%     firstnonan=startnonan(p)+1;
%     datnonan=datnonan(firstnonan:lastnonan);
% smoothamean=filtfilt(b,1,datnonan);
% amean=[nan(firstnonan-1,1); smoothamean'; nan(32-lastnonan,1)];
if ~isempty(varargin) && strcmp(varargin{1},'laminar')
    lineOut = plot(ax,amean_nan,F, 'Color', acolor,'linewidth',1.5); %% VERTICAL change color or linewidth to adjust mean line VERTICAL
else
    lineOut = plot(ax,F,amean_nan, 'Color', acolor,'linewidth',1.5); %% HORIZONTAL change color or linewidth to adjust mean line HORIZONTAL
end
if check
    hold off;
end

end


function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing
% need to average data before putting into code below
% fWidth = fWidth - 1 + mod(fWidth,2); %make sure filter length is odd
% dataStart = cumsum(dataIn(1:fWidth-2),2);
% dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
% dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
% dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
% dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
% dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];

dataOut=movmean(dataIn,fWidth,2);

end

