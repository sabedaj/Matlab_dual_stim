function [rate_all,Peak_latency_all]=spikelatencynoloops(chnnum,sp,IDs,trig,chn)
%determines number of spikes per ms 


TP = loadTrialParams;

rate_all=[];
Peak_latency_all=[];
    chnname=['Chn_' int2str(chnnum)];
        stimcchnname=['stimchn_' int2str(chn)];

        theseTrig = trig(cell2mat(TP(:,2)) == IDs(1,2))./30;
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
        rate_all.(chnname).(stimcchnname)(1:abs(BIN(1))+BIN(2)+1)=rate-mean(rate(1:80));
        [~,Peak_latency]=max(rate(abs(BIN(1))+1:end));
        Peak_latency=Peak_latency+abs(BIN(1));
        Peak_latency_all.(chnname).(stimcchnname)=Peak_latency;

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