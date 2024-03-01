function [V1_trialaverage,V2_trialaverage]=Ratesupression(ratestruct,savefilename)
%plot the average rate for the channels that are supressed and
%significantly surpressed

%iterate through stucture with folder entries and trial entries. 
% Calculate the average for each trial in each brain area. channels 1:64 is area V2. channel 65:128 is area V1.
% The structure of each trial is 128 channels x 181 time points x 40 repeats. 
% The average should be calculated for each channel and each time point. 
% Then the average should be calculated for each brain area and each trial within each folder.
%savefilename{1,1}{1,4} column 2 has the trial amplitude while column 1 has the trial number
uniquecurrents=unique(savefilename{1,1}{1,4}(:,2));%unique currents is a list of the unique currents used in the experiment
RS=cell(length(uniquecurrents),2);%RS stores averages for each unique current for the two areas
baseline=1:89;%baseline is the time point at which the average is calculated from
response_time=90+25:180;%response time is the time point at which the average is calculated from
for folder=1:length(ratestruct)
  
for trial=1:length(savefilename{folder}{1,4}(:,2))
    Trialnum=['T' num2str(trial)];
        %need to work out which channels are surpressed past 90ms compared to the first 90ms
        %then only include those channels in the average. Columnn 2 of ratestruct represents the time points
    %need to sort trials by amplitude and group them into a matrix for each amplitude
    currentOfTrial=savefilename{folder}{1,4}(trial,2);
    if all(currentOfTrial~=uniquecurrents)
        continue
    end

    temp1=mean(ratestruct{folder}.(Trialnum)(65:128,:,:),3);
    temp2=mean(ratestruct{folder}.(Trialnum)(1:64,:,:),3);


    Thresh_nobaseV1=mean(temp1(:,response_time),2)-(mean(temp1(:,baseline),2));%-1.*std(temp1(:,baseline),0,2));
    temp1(Thresh_nobaseV1>0,:)=NaN;%work out which are supressed
    Thresh_nobaseV2=mean(temp2(:,response_time),2)-(mean(temp2(:,baseline),2));%-1.*std(temp2(:,baseline),0,2));
    temp2(Thresh_nobaseV2>0,:)=NaN;%work out which are supressed
    RS{uniquecurrents==currentOfTrial,1}=[RS{uniquecurrents==currentOfTrial,1}; mean(temp1-(mean(temp1(:,baseline),2,'omitnan')),1,'omitnan').*1000];%V1 average
    RS{uniquecurrents==currentOfTrial,2}=[RS{uniquecurrents==currentOfTrial,2}; mean(temp2-(mean(temp2(:,baseline),2,'omitnan')),1,'omitnan').*1000];%V2 average

end
end
V1_trialaverage=RS(:,1);
V2_trialaverage=RS(:,2);


end