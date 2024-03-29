% EMBC 2023 analysis

% Heatmap alignment of stim electrode. Only include significant channels 

%sigmoids of V1 and V2

%latency plot V1 and V2

%no. resp/sig channels


%% loop through folders to get sigmoid data from each

D_data=dir;
namedat={D_data.name};
excludeforEMBC=[false false cellfun(@(x) (str2double(x(end-12:end-7))<221010),namedat(3:end),'UniformOutput',true)];%removing data for embc
D_data(excludeforEMBC,:)=[];%removing data for embc

sigmoidAll=cell(length(D_data),1);
rate_ALL=cell(length(D_data),1);
Peak_latencyALL=cell(length(D_data),1);
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
    [stimChn,~]=loadstimulationchannels;
    chn_range=1:64;

    if str2double(name(end-12:end-7))<220812
        warning('port D bad')
        chn_range(chn_range>96)=[];
    end
    timestart=2;
    timeend=8;% time to stop looking for spikes in ms
    [spk_all,rate_all,Peak_latency_all]=SigmoidGenerator(chn_range,sp,timestart,timeend);
    sigmoidAll{k}=spk_all;
    rate_ALL{k}=rate_all;
     Peak_latencyALL{k}=Peak_latency_all;
     %[rate,Peak_latency]=spikelatency(chn_range,sp);
%     spk_all_all={};
%     for chn=chn_range(1):chn_range(end)
%         chnname=['Chn_' num2str(chn)];
%         if mean(structfun(@mean, Peak_latency.(chnname)))-90<12
%             timestart=2;
%             timeend=22;% time to stop looking for spikes in ms
%         elseif mean(structfun(@mean, Peak_latency.(chnname)))-90>80
%             timestart=70;
%             timeend=90;% time to stop looking for spikes in ms
%         else
%             timestart=mean(structfun(@mean, Peak_latency.(chnname)))-10-90;%time after trig to start looking for spikes in ms
%             timeend=mean(structfun(@mean, Peak_latency.(chnname)))+10-90;% time to stop looking for spikes in ms
%         end
%         if ~isempty(sp{chn})
%             spk_all=SigmoidGenerator(chn,sp,timestart,timeend);
%         else
%             spk_all.(chnname)=nan(1,1);
%         end
%             
%         spk_all_all.(chnname)=spk_all.(chnname);
%     end
%     sigmoidAll{k}=spk_all_all;
%     rate_ALL{k}=rate;
%     Peak_latencyALL{k}=Peak_latency;
end
%% sigmoid
VA=2;%visual area recording
stimVA=1; %visual area being stimulated
cd([D_data(3).folder filesep D_data(3).name])
if VA==1
    columns=[5,7,8,6];
    chnnrange=65:128;
elseif VA==2
    columns=[1,3,4,2];
    chnnrange=1:64;
end
count=0;
sigmoidplot=nan(10,6);
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,columns);
if stimVA==2
stimchnarray=[1:16;33:48;49:64;17:32]';%+64;
elseif stimVA==1
    stimchnarray=[1:16;33:48;49:64;17:32]'+64;
end

for folder=1:length(sigmoidAll)
    currD = D_data(folder).name; % Get the current subdirectory name
    try
        cd([D_data(folder).folder filesep currD])
        if isempty(sigmoidAll{folder})
            continue;
        end
    catch
        continue
    end

    amp=loadAMP;
    for recchn=chnnrange(1):chnnrange(end)
        chnname=['Chn_' num2str(recchn)];
         [rrecchn,crecchn]=find(ordershapearray==recchn);
        fns = fieldnames(sigmoidAll{folder}.(chnname));
        for stimchn=1:length(fns)
            schn = str2double(regexp(fns{stimchn},'\d*','Match'));
            if schn>max(stimchnarray,[],'all') || schn<min(stimchnarray,[],'all') 
                continue
            end
             [r,c]=find(schn==stimchnarray);
%              if sigmoidAll{folder}.(chnname).(fns{stimchn})(end)>20
%                  load([currD(1:end-14) '.sp.mat'])
%                  ti=loadTrialInfo;
%                  ids=cell2mat(ti([false; cell2mat(ti(2:end,2))==sscanf(fns{stimchn},'stimchn_%d')],1));
%                  Online_raster(ids(end),sp{recchn})
%                  title([chnname ' ' fns{stimchn} ' ID_' num2str(ids(end))])
%              end
            %if r==rrecchn && c== crecchn
            %if sigmoidAll{folder}.(chnname).(fns{stimchn})(end)>0
                count=count+1;
                if size(sigmoidAll{folder}.(chnname).(fns{stimchn}),2)==length(amp)
                    sigmoidplot(count,amp)=sigmoidAll{folder}.(chnname).(fns{stimchn});
                else
                    sigmoidplot(count,[100; amp])=sigmoidAll{folder}.(chnname).(fns{stimchn});
                end
            %end
            %end
        end
    end
end
figure(3); hold on
if size(sigmoidplot,2)==100
errorbar([0; amp],nanmean(sigmoidplot(:,[100; amp])),nanstd(sigmoidplot(:,[100; amp]))./sqrt(sum(~isnan(sigmoidplot(:,end)))))
else
    errorbar([amp],nanmean(sigmoidplot(:,[amp])),nanstd(sigmoidplot(:,[amp]))./sqrt(sum(~isnan(sigmoidplot(:,end)))))
end
nanmean(sigmoidplot(:,[100; amp]))-max(nanmean(sigmoidplot(:,[100; amp])))/2;
xlabel('current (uA)')
ylabel('Firing rate (sp/s)')
set(gca,'TickDir','out');
axis square
%% rate
VA=2;%visual area recording
stimVA=1; %visual area being stimulated
cd([D_data(3).folder filesep D_data(3).name])
if VA==1
    columns=[5,7,8,6];
    chnnrange=65:128;
elseif VA==2
    columns=[1,3,4,2];
    chnnrange=1:64;
end
count=0;
rate_plot=nan(10,6);
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,columns);
if stimVA==2
stimchnarray=[1:16;33:48;49:64;17:32]';%+64;
elseif stimVA==1
    stimchnarray=[1:16;33:48;49:64;17:32]'+64;
end

for folder=1:length(rate_ALL)
    currD = D_data(folder).name; % Get the current subdirectory name
    try
        cd([D_data(folder).folder filesep currD])
        if isempty(rate_ALL{folder})
            continue;
        end
    catch
        continue
    end

    amp=loadAMP;
    for recchn=chnnrange(1):chnnrange(end)
        chnname=['Chn_' num2str(recchn)];
        [rrecchn,crecchn]=find(ordershapearray==recchn);
        fns = fieldnames(rate_ALL{folder}.(chnname));
        for stimchn=1:length(fns)
            schn = str2double(regexp(fns{stimchn},'\d*','Match'));
            if schn>max(stimchnarray,[],'all') || schn<min(stimchnarray,[],'all')
                continue
            end
            [r,c]=find(schn==stimchnarray);
            thresh=5;
            if any(rate_ALL{folder}.(chnname).(fns{stimchn})(end,2:8)>100)
                stop=0;
            end
%             if all(diff(sigmoidAll{folder}.(chnname).(fns{stimchn}))>=-10)
%                 ti=loadTrialInfo;
%                 tID=cell2mat(ti([false; cell2mat(ti(2:end,2))==schn],1));
%                 for i=1:length(tID)
%                 plotRawdata(recchn,tID(i))
%                 end
%             end
            %if mean(rate_ALL{folder}.(chnname).(fns{stimchn})(:,31+90))>thresh || mean(rate_ALL{folder}.(chnname).(fns{stimchn})(:,12+90))>thresh || mean(rate_ALL{folder}.(chnname).(fns{stimchn})(:,52+90))>thresh || mean(rate_ALL{folder}.(chnname).(fns{stimchn})(:,70+90))>thresh
               % continue%figure(50);hold on;plot(-90:90,mean(rate_ALL{folder}.(chnname).(fns{stimchn})))
           % end

            count=count+1;
                if size(rate_ALL{folder}.(chnname).(fns{stimchn}),1)==length(amp)
                    rate_plot(amp,1:181,count)=rate_ALL{folder}.(chnname).(fns{stimchn});
                else
                    rate_plot([100; amp],1:181,count)=rate_ALL{folder}.(chnname).(fns{stimchn});
                end

        end
    end
end

tempamp=[100; amp];
if VA==1
color1 = linspace(0,1,length(tempamp));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
else
    color1 = linspace(0,1,length(tempamp));
newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
end
figure; colororder(newcolors);
ax=axes;
hold on
F=-90:90;
for current=1:length(tempamp)
temp=squeeze(rate_plot(tempamp(current),:,:));
plottemp=nanmean(temp');
stdtemp=nanstd(temp')./sqrt(sum(~isnan(temp(1,:))));
fillOut = fill(ax,[F fliplr(F)],[plottemp+stdtemp fliplr(plottemp-stdtemp)],newcolors(current,:), 'FaceAlpha', 0.2,'linestyle','none');
plot(-90:90,[plottemp(1:89) nan nan nan  nan plottemp(94:181)])
end
legend('0','2','5','6','8','10')
set(gca,'TickDir','out');
%% sorting heatmap from sigmoiddata
VA=2;%visual area
stimVA=1; %visual area being stimulated
if VA==1
    columns=[5,7,8,6];
    chnnrange=65:128;
elseif VA==2
    columns=[1,3,4,2];
    chnnrange=1:64;
end
cd([D_data(3).folder filesep D_data(3).name])
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,columns);
if stimVA==2
stimchnarray=[1:16;33:48;49:64;17:32]';%+64;
elseif stimVA==1
    stimchnarray=[1:16;33:48;49:64;17:32]'+64;
end
heatmap_array_sigmoiddata=cell(100,1);
rate_distarray=cell(100,1);
%%%%%if the visual area is downstream or upstream
if VA~=stimVA
    sorted_array_sigmoiddata=cell(length(sigmoidAll),1);
    for folder=1:length(sigmoidAll)
        currD = D_data(folder).name; % Get the current subdirectory name
        try
            cd([D_data(folder).folder filesep currD])
            if isempty(sigmoidAll{folder})
                continue;
            end
        catch
            continue
        end
        for recchn=chnnrange(1):chnnrange(end)
            [rrecchn,crecchn]=find(ordershapearray==recchn);
            chnname=['Chn_' num2str(recchn)];
            fns = fieldnames(sigmoidAll{folder}.(chnname));
            for stimchn=1:length(fns)
                if size(sigmoidAll{folder}.(chnname).(fns{stimchn}),2)==length(amp)
                    for current=1:length(amp)
                        sorted_array_sigmoiddata{folder}{amp(current)}.(fns{stimchn})(rrecchn,crecchn)=sigmoidAll{folder}.(chnname).(fns{stimchn})(current);%iterates forwards
                    end
                else
                    ampmod=[100;amp];
                    for current=1:length(ampmod)
                        sorted_array_sigmoiddata{folder}{ampmod(current)}.(fns{stimchn})(rrecchn,crecchn)=sigmoidAll{folder}.(chnname).(fns{stimchn})(current);%iterates forwards
                    end
                end
            end
        end
    end
    
end
count=0;
count2=0;
heatmap_array=[];
rate_array=[];
for folder=1:length(sigmoidAll)
    currD = D_data(folder).name; % Get the current subdirectory name
    try
        cd([D_data(folder).folder filesep currD])
        if isempty(sigmoidAll{folder})
            continue;
        end
    catch
        continue
    end
    count=count+1;
for recchn=chnnrange(1):chnnrange(end)
    chnname=['Chn_' num2str(recchn)];
    [rrecchn,crecchn]=find(ordershapearray==recchn);
    
    fns = fieldnames(sigmoidAll{folder}.(chnname));
    for stimchn=1:length(fns)
        if recchn==chnnrange(1)
            for current=1:100
                heatmap_array_sigmoiddata{current}{count}.(fns{stimchn})=nan(31,7);
                rate_distarray{current}{count}.(fns{stimchn})=cell(31,7);
            end
        end
        schn=regexp(fns(stimchn),'\d*','Match');
        schn=str2double(schn{1});

        %if c~=crecchn %stim shank only
        if VA==stimVA
            [r,c]=find(schn==stimchnarray);
            if isempty(r)
                continue
            end
            rplus=16-r;
            cplus=4-c;
        end
            
            amp=loadAMP;
            if size(sigmoidAll{folder}.(chnname).(fns{stimchn}),2)==length(amp)
                for current=1:length(amp)
                    if VA~=stimVA
                        [r,c]=find(sorted_array_sigmoiddata{folder}{amp(current)}.(fns{stimchn})==max(sorted_array_sigmoiddata{folder}{amp(current)}.(fns{stimchn}),[],'all'));
                        rplus=16-r;
                        cplus=4-c;
                        if isempty(rplus)
                            rplus=1;
                            cplus=1;
                        end
                    end
                    heatmap_array_sigmoiddata{amp(current)}{count}.(fns{stimchn})(rrecchn+rplus(1),crecchn+cplus(1))=sigmoidAll{folder}.(chnname).(fns{stimchn})(current);%iterates forwards
                    rate_distarray{amp(current)}{count}.(fns{stimchn}){rrecchn+rplus(1),crecchn+cplus(1)}=rate_ALL{folder}.(chnname).(fns{stimchn})(current,:);
                end
            else
                tempamp=[100;amp];
                for current=1:length(tempamp)
                    if VA~=stimVA
                        [r,c]=find(sorted_array_sigmoiddata{folder}{tempamp(current)}.(fns{stimchn})==max(sorted_array_sigmoiddata{folder}{tempamp(current)}.(fns{stimchn}),[],'all'));
                        rplus=16-r;
                        cplus=4-c;
                        if isempty(rplus)
                            rplus=1;
                            cplus=1;
                        end
                    end
                    heatmap_array_sigmoiddata{tempamp(current)}{count}.(fns{stimchn})(rrecchn+rplus(1),crecchn+cplus(1))=sigmoidAll{folder}.(chnname).(fns{stimchn})(current);%iterates forwards
                    rate_distarray{tempamp(current)}{count}.(fns{stimchn}){rrecchn+rplus(1),crecchn+cplus(1)}=rate_ALL{folder}.(chnname).(fns{stimchn})(current,:);
                end
            end

        %end
        
        %         if sigmoidAll{folder}.(chnname).(fns{stimchn})(end-4)>100
        %             stop=0;
        %         end
        %         catch
        %             continue;
        %         end
        if recchn==chnnrange(end)
            count2=count2+1;
            for current=1:100
                heatmap_array{current}(:,:,count2)=heatmap_array_sigmoiddata{current}{count}.(fns{stimchn});
                rate_array{current}(:,:,count2)=rate_distarray{current}{count}.(fns{stimchn});
            end
        end
        
    end
end
end

%% DIstance, stimulation current effect
% current vs distance vs firing rate && distance vs firing rate

maxnum=1000;
step=100;
splitdist=cell(length(amp)+1,1);
splitratedist=splitdist;
splitdist(cellfun(@isempty, splitdist)) =   {nan((maxnum/step),1)};
stdsplitdist=splitdist;
countsplidist=splitdist;
splitratedist(cellfun(@isempty, splitratedist)) =   {nan((maxnum/step),181)};
stdratedist=splitratedist;
distcurrent=nan(20,6);
stddistcurrent=nan(20,6);
sigtest=zeros(maxnum/step,1);
sigtestRATE=zeros(91,1);

    color1 = linspace(0,1,length(tempamp));
    newcolors = [zeros(length(color1),1) color1' flipud(color1')];
    for current=1:length(tempamp)
        
        distanceArray=repmat(sqrt(([15:-1:0 1:15]'*50).^2+([3:-1:0 1:3]*200).^2),[1 1 size(heatmap_array{tempamp(current)},3)]);
        
        for i=step:step:maxnum-step
            
            splitdist{current}(i/step)=nanmean(heatmap_array{tempamp(current)}(distanceArray>=i-step & distanceArray<i));
            stdsplitdist{current}(i/step)=nanstd(heatmap_array{tempamp(current)}(distanceArray>=i-step & distanceArray<i))/sqrt(sum(~isnan(heatmap_array{tempamp(current)}(distanceArray>=i-50 & distanceArray<i))));
            countsplidist{current}(i/step)=nansum(~isnan(heatmap_array{tempamp(current)}(distanceArray>=i-step & distanceArray<i)));
            %             if current==1 %DISTANCE test significance between 0 and threshold value - 5uA for V2
            %                 sigtest(i/step)=ranksum(heatmap_array{tempamp(1)}(distanceArray>=i-step & distanceArray<i),heatmap_array{tempamp(2)}(distanceArray>=i-step & distanceArray<i));
            %             end
            splitratedist{current}(i/step,:)=mean(cell2mat(rate_array{tempamp(current)}(distanceArray>=i-step & distanceArray<i)));
            stdratedist{current}(i/step,:)=std(cell2mat(rate_array{tempamp(current)}(distanceArray>=i-step & distanceArray<i)))./sqrt(size(cell2mat(rate_array{tempamp(current)}(distanceArray>=i-step & distanceArray<i)),1));
            if current==1 %RATE test significance between 0 and threshold value - 5uA for V2
                current1rate=cell2mat(rate_array{tempamp(1)}(distanceArray>=0 & distanceArray<100));
                currentrate2=cell2mat(rate_array{tempamp(2)}(distanceArray>=0 & distanceArray<100));
                for timepoint=1:181
                    sigtestRATE(timepoint)=ranksum(current1rate(:,timepoint),currentrate2(:,timepoint));
                end
            end
            
            
        end
        
        figure(folder);colororder(newcolors);hold on; errorbar(0:step:maxnum-step,splitdist{current},stdsplitdist{current})
        xlabel('Approx. distance from the stim electrode(um)')
        ylabel('Firing rate (Sp/s)')
        %legend('2','5','6','8','10')
        
        %figure; plot(0:50:950,countsplidist); ylabel('# in average'); xlabel('Approx. distance from the stim electrode(um)')
    end
    set(gca,'TickDir','out');
    axis square
    legend('0','2','5','6','8','10')
    xlim([0, 800])
    figure(folder+50);hold on; plot(0:step:maxnum-step,countsplidist{1}); plot(0:step:maxnum-step,countsplidist{6})
  legend('0','2+')
  xlabel('Approx. distance from the stim electrode(um)')
  ylabel('# electrodes in average')
 axis square
% lightBLUE = [0 1 0];
% darkBLUE = [0 0 1];
% blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
% N=size(distcurrent,1);
%  figure
%  hold on
%  for dist=1:size(distcurrent,1)
%  errorbar([flipud(amp); 0], distcurrent(dist,:),stddistcurrent(dist,:), 'color',blueGRADIENTflexible(dist,N))
%  end
%  ylabel('Firing rate (Sp/s)')
%  xlabel('Current (ua)')
%  legendCell = cellstr(num2str([50:50:1000]', '%-d um'));
%  legend(legendCell)
amp_num=6;
if VA==1
color1 = linspace(0,1,size(splitratedist{amp_num},1)-1);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
else
    color1 = linspace(0,1,size(splitratedist{amp_num},1)-1);
newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
end

figure; colororder(newcolors);
ax=axes;
hold on
F=-90:90;
for i=1:size(splitratedist{amp_num},1)-1
    amean=splitratedist{amp_num}(i,:);
    astd=stdratedist{amp_num}(i,:);
fillOut = fill(ax,[F fliplr(F)],[amean+astd fliplr(amean-astd)],newcolors(i,:), 'FaceAlpha', 0.2,'linestyle','none');
end
fillOut.Annotation.LegendInformation.IconDisplayStyle = 'off';
colororder(newcolors);
plot(F,splitratedist{amp_num}(1:9,:));
legendCell = cellstr(num2str([0:step:maxnum-step]', '%-d'));
legend([legendCell;legendCell]);
xlabel('Time (ms)')
ylabel('Firing rate (sp/s)')
xline(0,'k')
xlim([-25,50])
set(gca,'TickDir','out');
axis square

%% average rate
for current=1:length(tempamp)
arrayratecat=nan(size(rate_array{tempamp(current)},1)*size(rate_array{tempamp(current)},2)*size(rate_array{tempamp(current)},3),181);
for i=1:size(rate_array{tempamp(current)},1)*size(rate_array{tempamp(current)},2)*size(rate_array{tempamp(current)},3)
    if ~isempty(rate_array{tempamp(current)}{i})
arrayratecat(i,:)=rate_array{tempamp(current)}{i};
    end

end
arrayratecat(isnan(arrayratecat(:,1)),:)=[];
figure(100); hold on; plot(-90:90,nanmean(arrayratecat,1))
end
   legend('0','2','5','6','8','10')
   xlabel('Time (ms)')
ylabel('Firing rate (sp/s)')
set(gca,'TickDir','out');
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
   try
    fns = fieldnames(sigmoidAll{folder}.(chnname));
   catch
       continue%no spikes found
   end
   
    for stimchn=1:length(fns)
        if recchn==chnnrange(1)
            for current=1:length(amp)+1
                heatmap_array_sigmoiddata{current}{count}.(fns{stimchn})=nan(16,4);
            end
        end

        
        if size(sigmoidAll{folder}.(chnname).(fns{stimchn}),2)==length(amp)
            for current=1:length(amp)
                heatmap_array_sigmoiddata{amp(current)}{count}.(fns{stimchn})(rrecchn,crecchn)=sigmoidAll{folder}.(chnname).(fns{stimchn})(current);%iterates forwards
            end
        else
            ampmod=[100;amp];
            for current=1:length(ampmod)
                heatmap_array_sigmoiddata{ampmod(current)}{count}.(fns{stimchn})(rrecchn,crecchn)=sigmoidAll{folder}.(chnname).(fns{stimchn})(current);%iterates forwards
            end
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

 amp_flip=[flipud(amp);0];
heatmap_centroid_centralised=cell(length(amp)+1,1);
 heatmap_centroid_centralised(cellfun(@isempty,heatmap_centroid_centralised)) = {nan(31,7)};
for current=1:length(amp)+1
    for penstim=1:size(heatmap_array{current},3)
        if all(isnan(heatmap_array{current}(:,:,penstim)),'all')
            continue
        end
        heatmap_it=heatmap_array{current}(:,:,penstim);
        heatmap_it(heatmap_it<0)=0;%make volume postive
        filtmov=ones(2,1)/2;

       filtered_heatmap=heatmap_it;%filtfilt(filtmov,1,heatmap_it);
        [r,c]=find(filtered_heatmap==max(filtered_heatmap,[],'all'));
        rplus=16-r;
        cplus=4-c;
        if isempty(rplus)
            rplus=1;
            cplus=1;
        end
        heatmap_centroid_centralised{current}(1+rplus:16+rplus,1+cplus:4+cplus,penstim)=filtered_heatmap;
    end

end
%%
  % DIstance, stimulation current effect
 % current vs distance vs firing rate && distance vs firing rate
splitdist=cell(length(amp)+1,1);
splitdist(cellfun(@isempty, splitdist)) =   {nan((1000/50),1)};
stdsplitdist=splitdist;
countsplidist=splitdist;
distcurrent=nan(20,6);
stddistcurrent=nan(20,6);
 for current=1:length(amp)+1

 distanceArray=repmat(sqrt(([15:-1:0 1:15]'*50).^2+([3:-1:0 1:3]*200).^2),[1 1 size(heatmap_centroid_centralised{current},3)]);

 for i=50:50:1000
     
    splitdist{current}(i/50)=nanmean(heatmap_centroid_centralised{current}(distanceArray>=i-50 & distanceArray<i));
     stdsplitdist{current}(i/50)=nanstd(heatmap_centroid_centralised{current}(distanceArray>=i-50 & distanceArray<i))/sum(~isnan(heatmap_centroid_centralised{current}(distanceArray>=i-50 & distanceArray<i)));
     countsplidist{current}(i/50)=nansum(~isnan(heatmap_centroid_centralised{current}(distanceArray>=i-50 & distanceArray<i)));
     distcurrent((i/50),current)=nanmean(heatmap_centroid_centralised{current}(distanceArray>=i-50 & distanceArray<i));
     stddistcurrent((i/50),current)=nanstd(heatmap_centroid_centralised{current}(distanceArray>=i-50 & distanceArray<i))/sum(~isnan(heatmap_centroid_centralised{current}(distanceArray>=i-50 & distanceArray<i)));
 end
 
 figure(1);hold on; errorbar(0:50:950,splitdist{current},stdsplitdist{current})
 xlabel('Approx. distance from the locus of neural activity(um)')
 ylabel('Firing rate (Sp/s)')
%legend('2','5','6','8','10')

%figure; plot(0:50:950,countsplidist); ylabel('# in average'); xlabel('Approx. distance from the stim electrode(um)')
 end
legend('10','8','6','5','2','0')
  %%
 %peak method
 amp_flip=[flipud(amp);0];
heatmap_centroid_centralised=cell(length(amp)+1,1);
 heatmap_centroid_centralised(cellfun(@isempty,heatmap_centroid_centralised)) = {nan(31,7)};
for current=1:length(amp)+1
    for penstim=1:size(heatmap_array{current},3)
        if all(isnan(heatmap_array{current}(:,:,penstim)),'all')
            continue
        end
        heatmap_it=heatmap_array{current}(:,:,penstim);
        heatmap_it(heatmap_it<0)=0;%make volume postive
        filtmov=ones(2,1)/2;

       filtered_heatmap=heatmap_it;%filtfilt(filtmov,1,heatmap_it);
        [r,c]=find(filtered_heatmap==max(filtered_heatmap,[],'all'));
        rplus=16-r;
        cplus=4-c;
        if isempty(rplus)
            rplus=1;
            cplus=1;
        end
        heatmap_centroid_centralised{current}(1+rplus:16+rplus,1+cplus:4+cplus,penstim)=filtered_heatmap;
    end
    figure; hHM=heatmap(nanmean(heatmap_centroid_centralised{current},3));
 xlabel('relative shank pos')
 ylabel('relative electrode pos')
 title(['Peak method ' num2str(amp_flip(current)) 'ua'])
end
 %% other methods of finding the locus of neural activity 
 % centroid method - nope
heatmap_centroid_centralised=cell(length(amp)+1,1);
 heatmap_centroid_centralised(cellfun(@isempty,heatmap_centroid_centralised)) = {nan(31,7)};
for current=1:length(amp)+1
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
    figure; hHM=heatmap(nanmean(heatmap_centroid_centralised{current},3));
 xlabel('relative shank pos')
 ylabel('relative electrode pos')
 title(['method 2:  centroid method ' num2str(amp_flip(current)) 'ua'])
end

 
 
%max area method - nope
heatmap_centroid_centralised=cell(length(amp)+1,1);
 heatmap_centroid_centralised(cellfun(@isempty,heatmap_centroid_centralised)) = {nan(31,7)};
for current=1:length(amp)+1
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
 
 
 %% linecut
 lightBLUE = [0 1 0];
darkBLUE = [0 0 1];
blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
avgFrmap=nanmean(heatmap_centroid_centralised{current},3);
std_heatmap=nanstd(heatmap_centroid_centralised{current},0,3)./sqrt(nansum(~isnan(heatmap_centroid_centralised{current}),3));
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
 
 figure
 plot([flipud(amp); 0],mean(distcurrent))
 xlabel('Current')
 ylabel('firing rate avg (sp/s)')
 
 %% rate and latency
 VA=2;%visual area
amp=loadAMP;

if VA==1
    columns=[5,7,8,6];
    chnnrange=65:128;
elseif VA==2
    columns=[1,3,4,2];
    chnnrange=1:64;
end
count=0;
rate_array=cell(length(sigmoidAll),1);
peak_array=cell(length(sigmoidAll),1);
 for folder=1:length(sigmoidAll)
     if isempty(sigmoidAll{folder})
         continue;
     end
    count=0;
     for recchn=chnnrange(1):chnnrange(end)
         chnname=['Chn_' num2str(recchn)];
         
         fns = fieldnames(sigmoidAll{folder}.(chnname));
         %count=count+1;
         for stimchn=1:length(fns)
             for i=1:size(rate_ALL{folder}.(chnname).(fns{stimchn}),1)
                 count=count+1;
                rate_array{folder}(count,:)=rate_ALL{folder}.(chnname).(fns{stimchn})(i,:);
                %peak_array{folder}{stimchn}(count)=Peak_latencyALL{folder}.(chnname).(fns{stimchn})(i);
%                 if ~isnan(rate_ALL{folder}.(chnname).(fns{stimchn})(i,1))
%                 plot(rate_ALL{folder}.(chnname).(fns{stimchn})(i,:))
%                 end
             end
         end
     end
 end
%   lightBLUE = [0 1 0];
% darkBLUE = [0 0 1];
% blueGRADIENTflexible = @(i,N) lightBLUE + (darkBLUE-lightBLUE)*((i-1)/(N-1));
% N=90;
%  figure; hold on
%  for stimchn=1:length(fns)
%      if ~isnan(nanmean(peak_array{folder}{stimchn}))
%  plot(-90:90,nanmean(rate_array{folder}{stimchn},1), 'color',blueGRADIENTflexible(nanmean(peak_array{folder}{stimchn})-90,N))
%      end
%  end
  for folder=1:length(sigmoidAll)
      if isempty(sigmoidAll{folder})
          continue;
      end
      figure; plot(-90:90,nanmean(rate_array{folder},1))
      title(num2str(folder))
  end
  %% heatmap of latency and peak value
hist_FR_lat=cell(10,1);
  for folder=1:length(sigmoidAll)-2
     
     if isempty(sigmoidAll{folder})
         continue;
     end
    count=0;
    fns = fieldnames(sigmoidAll{folder}.('Chn_1'));
    for stimchn=1:length(fns)
        hm_peaklat=nan(16,4);
        hm_FR=nan(16,4);
        for recchn=chnnrange(1):chnnrange(end)
            chnname=['Chn_' num2str(recchn)];
            [rrecchn,crecchn]=find(ordershapearray==recchn);
            
            count=count+1; 
            if ~isnan(Peak_latencyALL{folder}.(chnname).(fns{stimchn})(end))
            group=round((Peak_latencyALL{folder}.(chnname).(fns{stimchn})(end)-90)/10)+1;
            hist_FR_lat{group}=[hist_FR_lat{group} sigmoidAll{folder}.(chnname).(fns{stimchn})(end)];
            end
        end

    end
  end
  avgFRlat=zeros(9,1);
for groups=1:9
    avgFRlat(groups)=mean(hist_FR_lat{groups});
end
figure
bar(avgFRlat)