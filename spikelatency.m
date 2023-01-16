function [rate,Peak_latency]=spikelatency(spikeArray,ID)
%determines number of spikes per ms 
sp=spikeArray;
trig = loadTrig(0);
TP = loadTrialParams;
tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
theseTrig = trig(tID)./30;
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
 [~,Peak_latency]=max(rate);
 Peak_latency=Peak_latency+BIN(1);

end