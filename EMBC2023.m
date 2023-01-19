% EMBC 2023 analysis

% Heatmap alignment of stim electrode. Only include significant channels 

%sigmoids of V1 and V2

%latency plot V1 and V2

%no. resp/sig channels


%% loop through folders to get sigmoid data from each
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
%% sorting heatmap from sigmoiddata
VA=2;%visual area
amp=loadAMP;

if VA==1
    columns=[5,7,8,6];
    chnnrange=65:128;
elseif VA==2
    columns=[1,3,4,2];
    chnnrange=1:64;
end

order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,columns);
stimchnarray=[1:16;33:48;49:64;17:32]'+64;
heatmap_array_sigmoiddata=cell(length(amp)+1,1);
 %heatmap_array_sigmoiddata(cellfun(@isempty,heatmap_array_sigmoiddata)) = {nan(31,7)};
count=0;
count2=0;
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
            for current=1:length(amp)+1
                heatmap_array_sigmoiddata{current}{count}.(fns{stimchn})=nan(31,7);
            end
        end
        schn=regexp(fns(stimchn),'\d*','Match');
        schn=str2double(schn{1});
                [r,c]=find(schn==stimchnarray);
        rplus=16-r;
        cplus=4-c;

            for current=0:length(amp)
                if current==5 && length(sigmoidAll{folder}.(chnname).(fns{stimchn}))==5
                    continue
                end
                heatmap_array_sigmoiddata{current+1}{count}.(fns{stimchn})(rrecchn+rplus,crecchn+cplus)=sigmoidAll{folder}.(chnname).(fns{stimchn})(end-current);%iterates forwards
            end
            

%         if sigmoidAll{folder}.(chnname).(fns{stimchn})(end-4)>100
%             stop=0;
%         end
%         catch
%             continue;
%         end
        if recchn==chnnrange(end)
            count2=count2+1;
            for current=1:length(amp)+1
                heatmap_array{current}(:,:,count2)=heatmap_array_sigmoiddata{current}{count}.(fns{stimchn});
            end
        end
    
    end
end
end


%% heatmap plots
current=10;%ua
amp=loadAMP;
Indexcurrent=length(amp)-find(current==amp);%go from back since not all stimchn has zero included
 avgFrmap=nanmean(heatmap_array{Indexcurrent+1},3); figure; hHM=heatmap(avgFrmap);
 xlabel('relative shank pos')
 ylabel('relative electrode pos')
 title('Firing rate relative to stimulating electrode position (16,4)')
hHM.NodeChildren(3).YDir='normal'; 

lightBLUE = [0 1 0];
darkBLUE = [0 0 1];

blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
std_heatmap=nanstd(heatmap_array{Indexcurrent+1},0,3);%./nansum(~isnan(heatmap_array),3);
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
 
 g=heatmap_array{Indexcurrent+1};
 g(~isnan(g))=1;
 figure; heatmap(nansum(g, 3))
  xlabel('relative shank pos')
 ylabel('relative electrode pos')
title('#samples in average')


 
 %% DIstance, stimulation current effect
 % current vs distance vs firing rate && distance vs firing rate
splitdist=cell(length(amp)+1,1);
splitdist(cellfun(@isempty, splitdist)) =   {nan((1000/50),1)};
stdsplitdist=splitdist;
countsplidist=splitdist;
distcurrent=nan(20,6);
stddistcurrent=nan(20,6);
 for current=1:length(amp)+1

 distanceArray=repmat(sqrt(([15:-1:0 1:15]'*50).^2+([3:-1:0 1:3]*200).^2),[1 1 size(heatmap_array{current},3)]);

 for i=50:50:1000
     
    splitdist{current}(i/50)=nanmean(heatmap_array{current}(distanceArray>=i-50 & distanceArray<i));
     stdsplitdist{current}(i/50)=nanstd(heatmap_array{current}(distanceArray>=i-50 & distanceArray<i))/sum(~isnan(heatmap_array{current}(distanceArray>=i-50 & distanceArray<i)));
     countsplidist{current}(i/50)=nansum(~isnan(heatmap_array{current}(distanceArray>=i-50 & distanceArray<i)));
     distcurrent((i/50),current)=nanmean(heatmap_array{current}(distanceArray>=i-50 & distanceArray<i));
     stddistcurrent((i/50),current)=nanstd(heatmap_array{current}(distanceArray>=i-50 & distanceArray<i))/sum(~isnan(heatmap_array{current}(distanceArray>=i-50 & distanceArray<i)));
 end
 
 figure(18);hold on; errorbar(0:50:950,splitdist{current},stdsplitdist{current})
 xlabel('Approx. distance from the stim electrode(um)')
 ylabel('Firing rate (Sp/s)')
%legend('2','5','6','8','10')

%figure; plot(0:50:950,countsplidist); ylabel('# in average'); xlabel('Approx. distance from the stim electrode(um)')
 end
legend('10','8','6','5','2','0')
lightBLUE = [0 1 0];
darkBLUE = [0 0 1];
blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
N=size(distcurrent,1);
 figure
 hold on
 for dist=1:size(distcurrent,1)
 errorbar([flipud(amp); 0], distcurrent(dist,:),stddistcurrent(dist,:), 'color',blueGRADIENTflexible(dist,N))
 end
 ylabel('Firing rate (Sp/s)')
 xlabel('Current (ua)')
 legendCell = cellstr(num2str([50:50:1000]', '%-d um'));
 legend(legendCell)
 %% find v2 centroid

VA=2;%visual area record
VAstim=1;%visual are stim
amp=loadAMP;

if VA==1
    columns=[5,7,8,6];
    chnnrange=65:128;
elseif VA==2
    columns=[1,3,4,2];
    chnnrange=1:64;
end

order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,columns);
if VAstim==1
    stimchnarray=[1:16;33:48;49:64;17:32]'+64;
elseif VAstim==2
    stimchnarray=[1:16;33:48;49:64;17:32]';
end
heatmap_array_sigmoiddata=cell(length(amp)+1,1);
 %heatmap_array_sigmoiddata(cellfun(@isempty,heatmap_array_sigmoiddata)) = {nan(31,7)};
count=0;
count2=0;
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
            for current=1:length(amp)+1
                heatmap_array_sigmoiddata{current}{count}.(fns{stimchn})=nan(16,4);
            end
        end

            for current=0:length(amp)
                if current==5 && length(sigmoidAll{folder}.(chnname).(fns{stimchn}))==5
                    continue
                end
                heatmap_array_sigmoiddata{current+1}{count}.(fns{stimchn})(rrecchn,crecchn)=sigmoidAll{folder}.(chnname).(fns{stimchn})(end-current);%iterates forwards
            end
            
        if recchn==chnnrange(end)
            count2=count2+1;
            for current=1:length(amp)+1
                heatmap_array{current}(:,:,count2)=heatmap_array_sigmoiddata{current}{count}.(fns{stimchn});
            end
        end
    
    end
end
end

% centroid method
heatmap_centroid_centralised=cell(length(amp)+1,1);
 heatmap_centroid_centralised(cellfun(@isempty,heatmap_centroid_centralised)) = {nan(31,7)};
for current=1
    for penstim=1:size(heatmap_array{current},3)
        heatmap_it=heatmap_array{current}(:,:,penstim);
        heatmap_it(heatmap_it<0)=0;%make volume postive
        total_Volume=nansum(heatmap_it,'all');
        xcorr=nansum(heatmap_it,1);
         B=zeros(length(xcorr),1);
        for lims=2:length(xcorr)
            B(lims) =  trapz(1:lims, xcorr(1:lims));
        end
        [~,shankcentroid]=min(abs(B-(total_Volume/2)));
        
        ycorr=nansum(heatmap_it,2);
        C=zeros(length(ycorr),1);
        for lims=2:length(ycorr)
            C(lims) =  trapz(1:lims, ycorr(1:lims));
        end
        [~,electrodecentroid]=min(abs(C-(total_Volume/2)));
        rplus=16-electrodecentroid;
        cplus=4-shankcentroid;
        heatmap_centroid_centralised{current}(1+rplus:16+rplus,1+cplus:4+cplus,penstim)=heatmap_array{current}(:,:,penstim);
    end    
end
figure; hHM=heatmap(nanmean(heatmap_centroid_centralised{1},3));
 xlabel('relative shank pos')
 ylabel('relative electrode pos')
 title('method 2:  centroid method')
 
 
%max area method
heatmap_centroid_centralised=cell(length(amp)+1,1);
 heatmap_centroid_centralised(cellfun(@isempty,heatmap_centroid_centralised)) = {nan(31,7)};
for current=1
    for penstim=1:size(heatmap_array{current},3)
        heatmap_it=heatmap_array{current}(:,:,penstim);
        heatmap_it(heatmap_it<0)=0;%make volume postive
        total_Volume=nansum(heatmap_it,'all');
        xcorr=nansum(heatmap_it,1);
        [~,shankcentroid]=max(xcorr);
        
        ycorr=nansum(heatmap_it,2);
        [~,electrodecentroid]=max(ycorr);
        rplus=16-electrodecentroid;
        cplus=4-shankcentroid;
        heatmap_centroid_centralised{current}(1+rplus:16+rplus,1+cplus:4+cplus,penstim)=heatmap_array{current}(:,:,penstim);
    end    
end
figure; hHM=heatmap(nanmean(heatmap_centroid_centralised{1},3));
 xlabel('relative shank pos')
 ylabel('relative electrode pos')
 title('method 1:  max area method')
 
 %peak method
heatmap_centroid_centralised=cell(length(amp)+1,1);
 heatmap_centroid_centralised(cellfun(@isempty,heatmap_centroid_centralised)) = {nan(31,7)};
for current=1
    for penstim=1:size(heatmap_array{current},3)
        heatmap_it=heatmap_array{current}(:,:,penstim);
        heatmap_it(heatmap_it<0)=0;%make volume postive
       
        [r,c]=find(heatmap_it==max(heatmap_it,[],'all'));
        rplus=16-r;
        cplus=4-c;
        heatmap_centroid_centralised{current}(1+rplus:16+rplus,1+cplus:4+cplus,penstim)=heatmap_array{current}(:,:,penstim);
    end    
end
figure; hHM=heatmap(nanmean(heatmap_centroid_centralised{1},3));
 xlabel('relative shank pos')
 ylabel('relative electrode pos')
 title('method 1: peak method')
 
 
 
 
 %% linecut
 lightBLUE = [0 1 0];
darkBLUE = [0 0 1];
blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
avgFrmap=nanmean(heatmap_centroid_centralised{current},3);
std_heatmap=nanstd(heatmap_centroid_centralised{current},0,3)./nansum(~isnan(heatmap_centroid_centralised{current}),3);
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