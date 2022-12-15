function [samples_rand]=DownSample(pairavgcur,lengthneeded,s,seedpoint)

%%downsample
point=find(sum(~isnan(pairavgcur),2)>9,1,'first');% find the first row of pairs with 10+ columns on non-NAN values
firstsamples=find(~isnan(pairavgcur(point,:))); % select the first columns containing the 10+ non-NAN values
reset(s,seedpoint); %reset the seed to ensure consistency amongst all pairs of each current
if length(firstsamples)<lengthneeded %if there are less samples in firstsamples than we need to downsample
    samples_rand = datasample(s,1:size(pairavgcur,2),lengthneeded,'Replace',false); %randomly select enough samples to match lengthneeded out of the total number of pairs
    [~,pos]=intersect(samples_rand,firstsamples); %check how many random samples intersect with first samples
    samples_rand(:,pos)=[]; %remove intersecting points
    samples_rand=samples_rand(:,1:(lengthneeded-length(firstsamples))); %select the first random samples needed to increase firstsamples to lengthneeded
    samples_rand=[samples_rand firstsamples]; % concatenate the arrays
else
    datsamp = datasample(s,1:size(firstsamples,2),lengthneeded,'Replace',false);%if there are more samples in firstsamples, then we need to downsample and reduce it
    samples_rand=firstsamples(datsamp);
end
end