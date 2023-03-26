function sp=intersect_hist_tolerance(sp)
trig=ceil(loadTrig(0)./30+20000);%end recording
edges=0:1:trig(end);%edges of histogram bins
sp_count=zeros(128,trig(end));% holder for histogram spikes
for chn=1:128
[spx,~]=histcounts(sp{chn}(:,1),edges); %fit spikes into 1ms bins
sp_count(chn,:)=spx;
end
time_remove_spikes=find(sum(sp_count>0,1)>63);%find what time over 63 channels (i.e. one area) has spikes

for chn=1:128
    m=arrayfun(@(x)(x-1)<=sp{chn}(:,1) & sp{chn}(:,1)<x,time_remove_spikes,'UniformOutput',false);%find spike in original array to remove
    catm=horzcat(m{:});
    sp{chn}(sum(catm,2)>0,:)=[];%remove the spikes
end

end