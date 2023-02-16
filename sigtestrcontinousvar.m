function p=sigtestrcontinousvar(data,numresamp,varargin)
%numresamp= number of times to resample data
%data is the real data in columns
%varagin{1} determines tail, tail='both' is 2 tailed, 'left' or 'right' for
%single tail sig test. Default is 2 tailed
if nargin==3
    tail=varargin{1};
else
    tail='both';
end

numlabels=size(data,2);
numsamples=size(data,1);
%rng('default') % For reproducibility

slope=nan(numlabels,1);

for iteration=1:numresamp
        
    %randomise labels
    labelrandom=ceil(rand(numsamples,numlabels).*numlabels);
    
    %calculate slope
    fitdata = polyfit(labelrandom(:),data(:),1);
    slope(iteration)=fitdata(1);
        
end

%assumed gaussian distribution
%CI_low=mean(slope)-1.96*std(slope)/sqrt(numresamp);%confidence interval low
%CI_high=mean(slope)+1.96*std(slope)/sqrt(numresamp);%confidence interval high


labels=ones(size(data)).*(1:numlabels);%reat lavels for real data
fitdata = polyfit(labels(:),data(:),1);%fit real data

if (strcmp(tail,'both'))
    p1 = 1-mean(fitdata(1)<slope);
    p2 = 1-mean(fitdata(1)>slope);
    p=min(p1,p2);
elseif strcmp(tail,'left') 
    p = 1-mean(fitdata(1)<slope);
elseif strcmp(tail,'right')
    p = 1-mean(fitdata(1)>slope);
else
    p=1;
end
end