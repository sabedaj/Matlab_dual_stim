function p=SinglePairWErrorBars(ampinterest,plotYN)
% plot single pairs with trial repeats as error bars
ti=loadTrialInfo;
ti=cell2mat(ti(2:end,[1 2 18]));% trial ID, stimchn, AMP
E_MAP = Depth(6);
trialinterest=ti((ti(:,3)==ampinterest & ti(:,2)==0),1);% find trials

stimChn(1)=ti(ti(:,1)==trialinterest(1) & ti(:,2)~=0,2); % which stim chns
stimChn(2)=ti(ti(:,1)==trialinterest(2) & ti(:,2)~=0,2);

if stimChn(1)<17 %determines the shank with the stimulating electrodes
    shank=1;
elseif stimChn(1)<33 && stimChn(1)>16
    shank=4;%4 but for the purpose of removing stim shank==2
    stimChn=stimChn-16;
elseif stimChn(1)<49 && stimChn(1)>32
    shank=2;%2
    stimChn=stimChn-32;
else
    shank=3;%3
    stimChn=stimChn-48;
end

%load data
load('IDstruct.mat', 'IDstruct','baslinespikestruct')
if plotYN==1
    %setup figure
    vec = [100;80;50;30;15;0];
    hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
    N = 128;
    raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
    map = interp1(vec,raw,linspace(100,0,N),'pchip');
    cmap=colormap(map);
end
peaks=nan(40*4,5); %N_rep*shank and trial
for shankit=[1 3 4 2]
    if plotYN==1
        figure(shankit)
        hold on
        ax2=gca;
        hold(ax2, 'on');
        axes('Position',[ax2.Position(1) ax2.Position(2) ax2.Position(3) ax2.Position(4)*3/4])
        ax1=gca;
        ax2.Position=[ax2.Position(1) ax2.Position(2)+ax2.Position(4)*3/4 ax2.Position(3) ax2.Position(4)/4];
        hold(ax1,'off')
    end
    for trial=trialinterest(1):trialinterest(2)
        trialcheck=['T' num2str(trial)];
        baselinesubtracted_Data=IDstruct.(trialcheck)-baslinespikestruct.(trialcheck);
        baselinesubtracted_Data_sorted=baselinesubtracted_Data(E_MAP,:).*1000/(8-2);
        nantoadd=16-min(stimChn);
        dat2plot=[nan(nantoadd,size(baselinesubtracted_Data_sorted,2)); baselinesubtracted_Data_sorted(1+(16*(shankit-1)):16*shankit,:)];
        if plotYN==1
            [peak,~]=plotActivityAndPeak(dat2plot,abs(stimChn(1)-stimChn(2))-1,trial-trialinterest(1)+1,cmap,ax1,ax2,'across');
        else
            [peak,c]=find(dat2plot==max(dat2plot));
            peak(~diff(c))=[];
        end
        peaks(1+(size(baselinesubtracted_Data_sorted,2)*(shankit-1)):size(baselinesubtracted_Data_sorted,2)*shankit,trial-trialinterest(1)+1)=peak;
    end
    if plotYN==1
        xlim(ax1,[12 27])
        xlim(ax2,[12 27])
    end
end
if plotYN==0
    close all
end
peaks(sum(isnan(peaks),2)>0,:)=[];
p=friedman(peaks(:,2:4),1,'off');
%p=sigtestrcontinousvar(peaks(:,2:4),1000,'left');
end