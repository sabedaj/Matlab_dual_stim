% EMBC 2023 analysis

% Heatmap alignment of stim electrode. Only include significant channels 

%sigmoids of V1 and V2

%latency plot V1 and V2

%no. resp/sig channels


%%
D_data=dir;
totalSpikes4heatmap=cell(length(D_data),1);
totalSpikes4heatmap(cellfun(@isempty,totalSpikes4heatmap)) = {zeros(31,7)};
sigmoidAll=cell(length(D_data),1);
latency_peak_rate=cell(length(D_data),1);
latency_peak_rate(cellfun(@isempty,latency_peak_rate)) = {zeros(128,182)};
count4heatmap=cell(length(D_data),1);
count4heatmap(cellfun(@isempty,count4heatmap)) = {zeros(31,7)};
parfor k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    try
        cd([D_data(k).folder filesep currD])
        filepath=pwd;
        [filepathm,name,ext] = fileparts(filepath);
        %amp(k,1:5)=loadAMP;
        sp=load([name(1:end-14) '.sp.mat'],'sp');
        sp=sp.sp;
    catch
        continue
    end
    k
%     amplifier_channels=read_Intan_RHS2000_file;
%     ti=loadTrialInfo;
%     trial=cell2mat(ti([false; cell2mat(ti(2:end,18))==5],1));
    [stimChn,~]=loadstimulationchannels;
    if stimChn(1)<65
        VA=2;%visual area 2
        chn_range=1:64;
        columns=[1,3,4,2];
        stimchnarray=[1:16;33:48;49:64;17:32]';
    else
        VA=1;
        chn_range=1:128;
        columns=[1,3,4,2,5,7,8,6];
        stimchnarray=[1:16;33:48;49:64;17:32]'+64;
    end
    timestart=2;%time after trig to start looking for spikes in ms
    timeend=8;% time to stop looking for spikes in ms
%     order=Depth(1);
%     ordershapearray=reshape(order,16,8);
%     ordershapearray=ordershapearray(:,columns);
%     spikesArray=zeros(16,4);
    

   %%
%     for ID=trial'%5:5:size(ti,1)-1
%         schn=ti{ID,2}; %stimulaiton channel
%         [r,c]=find(schn==stimchnarray);
%         rplus=16-r;
%         cplus=4-c;
%         for chn=chn_range(1):chn_range(end)
%             %spikecount
%             NoSp=spikecount(timestart,timeend,ID,sp{chn}); %counts number of spieks after trig in time window
%             spikesArray(ordershapearray==chn)=NoSp;
%             [rrecchn,crecchn]=find(ordershapearray==chn);
%             totalSpikes4heatmap{k}(rrecchn+rplus,crecchn+cplus)=totalSpikes4heatmap{k}(rrecchn+rplus,crecchn+cplus)+NoSp;
%             count4heatmap{k}(rrecchn+rplus,crecchn+cplus)=count4heatmap{k}(rrecchn+rplus,crecchn+cplus)+1;
%             %latency
%             [rate,Peak_latency]=spikelatency(sp{chn},ID);
%             latency_peak_rate{k}(chn,:)=[Peak_latency,rate];
%             %sigmoid
%         end
%     end
    chanel_all=chn_range(1):chn_range(end);
    if str2double(name(end-12:end-7))<220812
        error('port D bad')
        chanel_all(chanel_all>96)=[];
    end
    spk_all=SigmoidGenerator(chanel_all,sp,timestart,timeend);
    sigmoidAll{k}=spk_all;
end
%%
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
%% plotting heatmap from sigmoiddata
order=Depth(1);
ordershapearray=reshape(order,16,8);
columns=[5,7,8,6];
ordershapearray=ordershapearray(:,columns);
stimchnarray=[1:16;33:48;49:64;17:32]'+64;
 heatmap_array_sigmoiddata=cell(24,1);
 %heatmap_array_sigmoiddata(cellfun(@isempty,heatmap_array_sigmoiddata)) = {nan(31,7)};
count=0;
count2=0;
chnnrange=65:128;
heatmap_array=[];
for folder=1:length(sigmoidAll)
    if isempty(sigmoidAll{folder})
        continue;
    end
    count=count+1;
for recchn=chnnrange(1):chnnrange(end)
    chnname=['Chn_' num2str(recchn)];
    [rrecchn,crecchn]=find(ordershapearray==recchn);
    
    fns = fieldnames(sigmoidAll{folder}.(chnname));
    for stimchn=1:length(fns)
        if recchn==chnnrange(1)
            heatmap_array_sigmoiddata{count}.(fns{stimchn})=nan(31,7);
        end
        schn=regexp(fns(stimchn),'\d*','Match');
        schn=str2double(schn{1});
                [r,c]=find(schn==stimchnarray);
        rplus=16-r;
        cplus=4-c;
        try
        heatmap_array_sigmoiddata{count}.(fns{stimchn})(rrecchn+rplus,crecchn+cplus)=sigmoidAll{folder}.(chnname).(fns{stimchn})(end-4);
%         if sigmoidAll{folder}.(chnname).(fns{stimchn})(end-4)>100
%             stop=0;
%         end
        catch
            continue;
        end
        if recchn==chnnrange(end)
            count2=count2+1;
            heatmap_array(:,:,count2)=heatmap_array_sigmoiddata{count}.(fns{stimchn});
        end
    
    end
end
end
%%
 avgFrmap=nanmean(heatmap_array,3); figure; hHM=heatmap(avgFrmap);
 xlabel('relative shank pos')
 ylabel('relative electrode pos')
 title('Firing rate relative to stimulating electrode position (16,4)')
hHM.NodeChildren(3).YDir='normal'; 

lightBLUE = [0 1 0];
darkBLUE = [0 0 1];

blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
std_heatmap=nanstd(heatmap_array,0,3);%./nansum(~isnan(heatmap_array),3);
figure;
hold on
N=7;
for i=1:7
%plot(avgFrmap(:,i),'color',blueGRADIENTflexible(i,N))
errorbar(1:31,avgFrmap(:,i),std_heatmap(:,i),'color',blueGRADIENTflexible(i,N))
end
 xlabel('relative electrode')
 ylabel('Firing rate (sp/s)')
 legend
 
 g=heatmap_array;
 g(~isnan(g))=1;
 figure; heatmap(nansum(g, 3))
  xlabel('relative shank pos')
 ylabel('relative electrode pos')
 
title('#samples in average')
 
 %%
 distanceArray=repmat(sqrt(([15:-1:0 1:15]'*50).^2+([3:-1:0 1:3]*200).^2),[1 1 size(heatmap_array,3)]);
 splitdist=nan(1000/50,1);
stdsplitdist=nan(1000/50,1);
countsplidist=nan(1000/50,1);
 for i=50:50:1000
     
    splitdist(i/50)=nanmean(heatmap_array(distanceArray>=i-50 & distanceArray<i));
     stdsplitdist(i/50)=nanstd(heatmap_array(distanceArray>=i-50 & distanceArray<i))/sum(~isnan(heatmap_array(distanceArray>=i-50 & distanceArray<i)));
     countsplidist(i/50)=nansum(~isnan(heatmap_array(distanceArray>=i-50 & distanceArray<i)));
 end
 figure(18);hold on; errorbar(0:50:950,splitdist,stdsplitdist)
 xlabel('Approx. distance from the stim electrode(um)')
 ylabel('Firing rate (Sp/s)')
%legend('2','5','6','8','10')
legend('10','8','6','5','2')
%figure; plot(0:50:950,countsplidist); ylabel('# in average'); xlabel('Approx. distance from the stim electrode(um)')


%% find v2 centroid
