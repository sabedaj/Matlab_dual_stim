function corrV1V2=spike_count_corr(spksV1,spksV2,V1time,V2time)
%spksV1 is a 64 channelx181 timepoint x N Trials matrix of the spike count at each timepoint for each channel
%spksV2 is a N trial x 181 timepoint matrix of the spike count at each timepoint 
%collapse the V1 array into N trials x 181 timepoints
reshaped_dataV1 = reshape(spksV1, 64, [], size(spksV1, 2));
meanV1 = squeeze(mean(reshaped_dataV1, 1,'omitnan'));
meanV1=meanV1-mean(meanV1(:,1:89),2,'omitnan');
spksV2=spksV2-mean(spksV2(:,1:89),2,'omitnan');
%calculate the correlation between the spike count from 92:112ms in V1 and 92:180ms in V2
corrV1V2 = corr(sum(meanV1(:,V1time),2),sum(spksV2(:,V2time),2));

end