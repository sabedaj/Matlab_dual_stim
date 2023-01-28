function [rate_all,Peak_latency_all]=spikelatency(chnnum,spikeArray,varargin)
%determines number of spikes per ms 


TP = loadTrialParams;
if nargin==2
    trialinfo=loadTrialInfo;
    trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
    numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
    chn=unique(trialinfo(:,2));
    IDs=cell(length(chn),1);
    trig = loadTrig(0);
    count=0;
    for trialloop=1:numsim_elect:length(trialinfo(:,2))
        if sum(trialinfo(trialloop:trialloop+numsim_elect-1,2)==0)==numsim_elect-1
            count=count+1;
            for chnloop=1:length(chn)
                if any(trialinfo(trialloop:trialloop+numsim_elect-1,2)==chn(chnloop))
                    IDs{chnloop}(count,1:2)=[trialinfo(trialloop,17) trialloop];
                end
            end
            %allsingletrialIDs(count)=trialloop;%index   trialinfo(trialloop,1);
        end
    end
else

    trig=varargin{2};
    chn=varargin{3};
    IDs{1}=varargin{1};%not this is only if you want to put a single ID in with a single stimchn
end
rate_all=[];
Peak_latency_all=[];
for recordchn=1:length(chnnum)
    chnname=['Chn_' int2str(chnnum(recordchn))];
    sp=spikeArray{chnnum(recordchn)};
    for stimchn=1:length(chn)
        stimcchnname=['stimchn_' int2str(chn(stimchn))];
        if nargin==2
            IDs{stimchn}=sortrows(IDs{stimchn});
            IDs{stimchn}=IDs{stimchn}(any(IDs{stimchn},2),:);
        end

        for ID=1:size(IDs{stimchn},1)
            %t=['ID_' num2str((IDs{stimchn}(ID,2)))];
            tID = find(cell2mat(TP(:,2)) == IDs{stimchn}(ID,2));
            theseTrig = trig(tID)./30;
            theseTrig(theseTrig<0)=[];
            nT=length(theseTrig);
            BIN = [-90 90];
            xdata = [];
            
            for tr = 1:nT
                theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
                for i = 1:length(theseSp)
                    xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
                end
            end
            Z = hist(xdata,0:length(BIN(1):BIN(2)-1)); %#ok<HIST>
            rate = (1000/nT)*Z;
            rate_all.(chnname).(stimcchnname)(ID,1:abs(BIN(1))+BIN(2)+1)=rate-mean(rate(1:80));
            [~,Peak_latency]=max(rate(abs(BIN(1))+1:end));
            Peak_latency=Peak_latency+abs(BIN(1));
            Peak_latency_all.(chnname).(stimcchnname)(ID)=Peak_latency;
        end
    end
end
% 
% for chn=chnnum
% sp=spikeArray{chn};
% trig = loadTrig(0);
% TP = loadTrialParams;
% tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
% theseTrig = trig(tID)./30;
% nT=length(theseTrig);
% 
% BIN = [-90 90]; 
% xdata = [];
% 
% for tr = 1:nT
%     theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
%     for i = 1:length(theseSp)
%         xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
%     end
% end
%  Z = hist(xdata,0:length(BIN(1):BIN(2)-1)); %#ok<HIST>
%  rate = (1000/nT)*Z;
%  [~,Peak_latency]=max(rate);
%  Peak_latency=Peak_latency+BIN(1);
% 
% end
end