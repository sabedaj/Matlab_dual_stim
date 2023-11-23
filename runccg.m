timeinms=300000;
window=50;
ccgall=zeros(1,window*2);
count=0;
for chn1=65:128
spikeTrain1=sp{chn1}(sp{chn1}(:,1)<timeinms,1);
if length(spikeTrain1)>150
    figure
    plot(sp{chn1}(sp{chn1}(:,1)<timeinms,2:end)')
end
for chn2= 1:64
spikeTrain2=sp{chn2}(sp{chn2}(:,1)<timeinms,1);
[ccg, lags] = crossCorrelogram(spikeTrain2,spikeTrain1,window);

%check spikes are good
if length(spikeTrain2)>150 && chn1==1
    figure
    plot(sp{chn2}(sp{chn2}(:,1)<timeinms,2:end)')
end

ccg = ccg / ((length(spikeTrain1)+length(spikeTrain2)));%normalize by number of spikes in train1
if ~isempty(ccg)&& any(ccg>mean(ccg)+std(ccg)*4)
ccgall=ccgall+ccg;
count=count+1;
end
if any(ccg>0.25) && ~isempty(ccg) && (length(spikeTrain1)+length(spikeTrain2)>10)
figure
bar(lags(1:end-1), ccg);
title([num2str(chn1) ' & ' num2str(chn2)])
end
end
end
figure
bar(lags(1:end-1), ccgall*100./count);
title('Basline CCG between V1 and V2 summed over all channels and normalised')
xlabel('Time(ms)')
ylabel('Percentage of V1 spikes correlated with V2 spikes')
set(gca,'TickDir','out');