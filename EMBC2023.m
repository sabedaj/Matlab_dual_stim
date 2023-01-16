% EMBC 2023 analysis

% Heatmap alignment of stim electrode. Only include significant channels 

%sigmoids of V1 and V2

%latency plot V1 and V2

%no. resp/sig channels

D_data=dir;
totalSpikes4heatmap=cell(length(D_data),1);
totalSpikes4heatmap(cellfun(@isempty,totalSpikes4heatmap)) = {zeros(31,7)};
count4heatmap=cell(length(D_data),1);
count4heatmap(cellfun(@isempty,count4heatmap)) = {zeros(31,7)};
parfor k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    try
        cd([D_data(k).folder filesep currD])
        filepath=pwd;
        [filepathm,name,ext] = fileparts(filepath);
        sp=load([name(1:end-14) '.sp.mat'],'sp');
    catch
        continue
    end
    k
    amplifier_channels=read_Intan_RHS2000_file;
    ti=loadTrialInfo;
    trial=cell2mat(ti([false; cell2mat(ti(2:end,18))==5],1));
    [stimChn,~]=loadstimulationchannels;
    if stimChn(1)<65
        VA=2;%visual area 2
        chn_range=1:64;
        columns=[1,3,4,2];
        stimchnarray=[1:16;33:48;49:64;17:32]';
    else
        VA=1;
        chn_range=65:128;
        columns=[5,7,8,6];
        stimchnarray=[1:16;33:48;49:64;17:32]'+64;
    end
    timestart=2;%time after trig to start looking for spikes in ms
    timeend=8;% time to stop looking for spikes in ms
    order=Depth(1);
    ordershapearray=reshape(order,16,8);
    ordershapearray=ordershapearray(:,columns);
    spikesArray=zeros(16,4);
    

   %%
    for ID=trial'%5:5:size(ti,1)-1
        schn=ti{ID,2}; %stimulaiton channel
        [r,c]=find(schn==stimchnarray);
        rplus=16-r;
        cplus=4-c;
        for chn=chn_range(1):chn_range(end)
            %spikecount
            NoSp=spikecount(timestart,timeend,ID,sp.sp{chn}); %counts number of spieks after trig in time window
            spikesArray(ordershapearray==chn)=NoSp;
            [rrecchn,crecchn]=find(ordershapearray==chn);
            totalSpikes4heatmap{k}(rrecchn+rplus,crecchn+cplus)=totalSpikes4heatmap{k}(rrecchn+rplus,crecchn+cplus)+NoSp;
            count4heatmap{k}(rrecchn+rplus,crecchn+cplus)=count4heatmap{k}(rrecchn+rplus,crecchn+cplus)+1;
            %latency
            [rate,Peak_latency]=spikelatency(spikeArray,ID);
            %sigmoid
        end
    end

end

counter=0;
x=zeros(31,7,3);
for i=1:length(count4heatmap)
    if sum(totalSpikes4heatmap{i},'all')~=0
        counter=counter+1;
        x(:,:,counter)=totalSpikes4heatmap{i}./count4heatmap{i};
    end
end
y=nanmean(x,3);
figure; heatmap(y)