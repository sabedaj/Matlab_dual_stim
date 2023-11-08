window=200;
ccgall=zeros(1,window*2);
for chn1=65:128
spikeTrain1=sp{chn1}(sp{chn1}(:,1)<20000,1);
if length(spikeTrain1)>150
    figure
    plot(sp{chn1}(sp{chn1}(:,1)<20000,2:end)')
end
for chn2= 1:64
spikeTrain2=sp{chn2}(sp{chn2}(:,1)<20000,1);
[ccg, lags] = crossCorrelogram(spikeTrain1,spikeTrain2,window);
if length(spikeTrain2)>150 && chn1==1
    figure
    plot(sp{chn2}(sp{chn2}(:,1)<20000,2:end)')
end
ccg = ccg / ((length(spikeTrain1)+length(spikeTrain2)));%normalize by number of spikes in train1
ccgall=ccgall+ccg;
if any(ccg>0.45) && ~isempty(ccg) && (length(spikeTrain1)+length(spikeTrain2)>10)
figure
bar(lags(1:end-1), ccg);
title([num2str(chn1) ' & ' num2str(chn2)])
end
end
end
figure
bar(lags(1:end-1), ccgall);