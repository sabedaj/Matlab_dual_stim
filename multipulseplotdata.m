function multipulseplotdata(ratestruct,savefilename,area)
s = RandStream('mt19937ar','Seed',296);
AMPall=[2 5 6 8 10];
 map2=[255,207,0;255,140,143;219,125,255;140,56,255;13,0,102]/255;
 map1=[181,2,0;255,140,74;255,207,0;77,217,255;36,77,255]/255;
 map4=[255,176,0;254,97,0;220,38,127;112,66,255;77,199,255]/255;
 map3=[128,0,255;195,164,207;120,120,120;18,125,9;128,255,0]/255;
baselinetime=1:89;
Pulseall=1:5;
 subselectbaseline=1:2:89;
 SD_supression=0.5;
 % cutoffs for each pen
 startpen=[datetime(2022,09,06,22,20,00) datetime(2022,09,07,19,30,00) datetime(2022,10,10,18,30,00) datetime(2022,10,11,16,30,00) datetime(2022,11,21,18,30,00) datetime(2022,11,22,22,00,00)];
 endpen=[datetime(2022,09,07,18,30,00) datetime(2022,09,08,15,30,00) datetime(2022,10,11,13,30,00) datetime(2022,10,12,11,30,00) datetime(2022,11,22,16,30,00) datetime(2022,11,23,12,00,00)];
 
% pre-processing
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=[];
ratespiking_pensplit=[];
Ratestitchspiking=[];
windowsize=3.3333;%in ms for stitching the data 3.33333
if strcmp(area,'v1')
    chnrng=65:128;
elseif strcmp(area,'v2')
    chnrng=1:64;
end

for numerfolders=1:numfolderstotal
    trialend=size(savefilename{numerfolders}{4},1);
    yrmnth=['YM' savefilename{numerfolders}{3}(1:6)];
  %need to loop through and find the maximum response of any trial for each channel
  maxdatchn=zeros(64,1);
  for trial=1:trialend
      trialcheck=['T', num2str(trial)];
    dat=max(nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,1:181,:).*1000,3)-mean(nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,1:90,:).*1000,3),2));
    for chn=1:64
        if dat(chn)>maxdatchn(chn)
            maxdatchn(chn)=dat(chn);
        end
    end
  end


    for trial=1:trialend
        trialcheck=['T', num2str(trial)];
        stimchncheck=['SChn', num2str(savefilename{numerfolders}{4}(trial,3))];
        AMPcheck=['A' num2str(savefilename{numerfolders}{4}(trial,2))];
        if strcmp(AMPcheck,'A-1')
            continue; %AMPcheck='A0';
        end
        Pulse=savefilename{numerfolders}{4}(trial,4);
        Pulsecheck=['P' num2str(Pulse)];
        fold_int=dir(['*' savefilename{numerfolders}{5}]);
        %C:\data\multipulse\PEN3_V1stim_221122_222759 chn=28+64, trial 14
        %no longer the one below
%         if strcmp(savefilename{numerfolders}{5},'220907_235749') && trial==15
%             %make first plot showing individual responses to stim
%             cd(fold_int.name);
%             plotstackedtrigs(70, trial, (1:40), 'DT')
%             cd(fold_int.folder);
%         end
        %stitching
        
        chnstitch=['Chn' num2str(savefilename{numerfolders}{4}(trial,3))];

        for pulsit=5:-1:Pulse+1
            Pulsecheckstitch=['P' num2str(pulsit)];
            if pulsit>1
            Ratestitchspiking.(yrmnth).(AMPcheck).(chnstitch).(Pulsecheckstitch)(:,94+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse))=nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,94+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse),:).*1000,3);
            else
                Ratestitchspiking.(yrmnth).(AMPcheck).(chnstitch).(Pulsecheckstitch)(:,91+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse))=nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,91+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse),:).*1000,3);
            end
        end
        if Pulse>1
        Ratestitchspiking.(yrmnth).(AMPcheck).(chnstitch).(Pulsecheck)(:,94+ceil(windowsize*(Pulse-1)):391)=nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,94+ceil(windowsize*(Pulse-1)):391,:).*1000,3);
        else
            Ratestitchspiking.(yrmnth).(AMPcheck).(chnstitch).(Pulsecheck)(:,91+ceil(windowsize*(Pulse-1)):391)=nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,91+ceil(windowsize*(Pulse-1)):391,:).*1000,3);
        end
        datatoinclude=nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,:,:).*1000,3)-mean(nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,1:90,:).*1000,3),2);
          %datatoinclude=datatoinclude./maxdatchn;
    
        %datatoinclude=datatoinclude./max(datatoinclude(:,1:181),[],2);%-nanmean(datatoinclude(:,1:90),2))./std(datatoinclude(:,1:90),[],2);%z-score
        datatoinclude(isinf(datatoinclude))=nan;
        if isfield(ratespiking,(AMPcheck)) && isfield(ratespiking.(AMPcheck),(Pulsecheck))
            ratespiking.(AMPcheck).(Pulsecheck)=cat(1,ratespiking.(AMPcheck).(Pulsecheck),datatoinclude);%-nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,baselinetime,:).*1000,[2 3]));
            flnm.(AMPcheck).(Pulsecheck)=[flnm.(AMPcheck).(Pulsecheck);savefilename{numerfolders}{5}];%-nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,baselinetime,:).*1000,[2 3]));
        else
            ratespiking.(AMPcheck).(Pulsecheck)=datatoinclude;
            flnm.(AMPcheck).(Pulsecheck)=[savefilename{numerfolders}{5}];
        end
        %data grouped based on penetration
        % first check whether the datetime is within pen range  savefilename{numerfolders}{5} is in yymmdd_hhmmss format ensure pen is at the front 
        for pen=1:length(startpen)
            if datetime(savefilename{numerfolders}{5},'InputFormat','yyMMdd_HHmmss')>startpen(pen) && datetime(savefilename{numerfolders}{5},'InputFormat','yyMMdd_HHmmss')<endpen(pen)
                pencheck=['Pen' num2str(pen)];
                if pen==4 && ~isempty(fold_int) && ~strcmp(fold_int.name(1:10),'multipulse')
                    continue
                end
                if isfield(ratespiking_pensplit,(pencheck)) && isfield(ratespiking_pensplit.(pencheck),(AMPcheck)) && isfield(ratespiking_pensplit.(pencheck).(AMPcheck),(stimchncheck))&& isfield(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck),(Pulsecheck))
                    ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)=cat(1,ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck),datatoinclude);%-nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,baselinetime,:).*1000,[2 3]));
                else
                    ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)=datatoinclude;
                end
            end
        end

        
         thresh=squeeze(mean(datatoinclude(:,subselectbaseline),2)+(std(datatoinclude(:,subselectbaseline),[],2).*3));
        notsigsingle=squeeze(mean(datatoinclude(:,92:181),2))<=thresh;
        Ratestitchspiking.(yrmnth).(AMPcheck).(chnstitch).(Pulsecheck)(notsigsingle,:)=nan;

    end
end
ratespiking_suppression=ratespiking;%use this later to check suppression


%plot significant values for all recorded channels
ratespiking_mean=zeros(length(AMPall),length(Pulseall));
ratespiking_std=zeros(length(AMPall),length(Pulseall));
ratepeak=zeros(length(AMPall),length(Pulseall));
ratepeakstd=zeros(length(AMPall),length(Pulseall));
numsig=zeros(length(AMPall),5);
numelect=zeros(length(AMPall)+1,5+1);
numelect(1,2:end)=Pulseall;
numelect(2:end,1)=AMPall';
StimElectcount=zeros(length(AMPall),5);
%loop through amp and number of pulses then plot
avgtime=91:130;
figure
	p = panel();
	p.pack(2, 3);
for AMP=1:length(AMPall)
    AMPcheck=['A' num2str(AMPall(AMP))];
    if AMP<4
        p(1,AMP).select();
    else
        p(2,AMP-3).select();
    end
    ax=p.de.axis;
    ax=ax(AMP);
    hold on
    
    for Pulse=1:length(Pulseall)
        Pulsecheck=['P' num2str(Pulseall(Pulse))];
        %make sure only responding elect included
       timepeak=round(90+Pulse*3.333):round(120+Pulse*3.333);
        if isfield(ratespiking.(AMPcheck),(Pulsecheck))
        thresh=squeeze(mean(ratespiking.(AMPcheck).(Pulsecheck)(:,subselectbaseline),2)+(std(ratespiking.(AMPcheck).(Pulsecheck)(:,subselectbaseline),[],2).*3));
        notsigsingle=squeeze(mean(ratespiking.(AMPcheck).(Pulsecheck)(:,timepeak),2))<=0;
        ratespiking.(AMPcheck).(Pulsecheck)(notsigsingle,:)=nan;
        numsig(AMP,Pulse)=sum(~notsigsingle);
        %ratespiking.(AMPcheck).(Pulsecheck)= ratespiking.(AMPcheck).(Pulsecheck)-squeeze(mean(ratespiking.(AMPcheck).(Pulsecheck)(:,subselectbaseline),2));
        dat=ratespiking.(AMPcheck).(Pulsecheck);
        numelect(AMP+1,Pulse+1)= sum(~isnan(dat(:,1)));
        % Reshape the array into a 2D matrix with 128 elements per column
        reshapedArray = reshape(dat(:,1), 64, []);
        % Count the number of non-NaN values in each column
        StimElectcount(AMP,Pulse) = sum(~isnan(sum(~isnan(reshapedArray))));
        %dat=dat./max(dat(:,1:181),[],2);
        maxdat=max(dat(:,round(90+Pulse*3.333):round(95+Pulse*3.333)),[],2);
        
        ratespiking_mean(AMP,Pulse)=nanmean(dat(:,timepeak),'all');%mean(sum(ratespiking.(AMPcheck).(Pulsecheck)(:,timepeak),2),'omitnan');%%nanmean(maxdat(maxdat~=0),'all');%nanmean(ratespiking.(AMPcheck).(Pulsecheck)(:,avgtime),'all');
        ratespiking_std(AMP,Pulse)=SEM(dat(:,timepeak),0);
        %peakvals
        ratepeak(AMP,Pulse)=mean(maxdat,'omitnan');
        ratepeakstd(AMP,Pulse)=SEM(maxdat,0);
       
        stdshade(dat,0.2,[0.0784+(Pulse-1)*(1-0.0784)/(length(Pulseall)-1) 0.5647+(Pulse-1)*(0.8431-0.5647)/(length(Pulseall)-1) 1-(Pulse-1)*(1-0)/(length(Pulseall)-1)],[-90:300],1,ax)
        end
    end
    lgd=legend('1','2','3','4','5');
    lgd.Title.String = '# pulses';
    xlim([-90,90])
    axis square
    title(['All data' num2str(AMPall(AMP)) '\muA'])
    ylabel('Sp/s')
    xlabel('Time (ms)')
    set(gca,'TickDir','out');
    ylim([-0.5 1])
    for j = 2:size(numelect,2)
        text(j*10-80, 65, num2str(numelect(AMP+1,j)), 'HorizontalAlignment', 'center');
    end
end



% p(2,3).select();
% hold on;
% for Pulse=1:length(Pulseall)
% errorbar(AMPall,ratespiking_mean(:,Pulse),ratespiking_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
% end
%  lgd=legend('1','2','3','4','5');
% lgd.Title.String = '# pulses';
%     title('Average 1-40ms after stim')
%     ylabel('Sp/s')
%     xlabel('Current \muA')
%     set(gca,'TickDir','out');
    
    figure;
    p = panel;
    p.pack(1,2)
    p(1,1).select();
    hold on
    for Pulse=1:length(Pulseall)
        errorbar(AMPall,ratepeak(:,Pulse),ratepeakstd(:,Pulse),'color',[0.0784+(Pulse-1)*(1-0.0784)/(length(Pulseall)-1) 0.5647+(Pulse-1)*(0.8431-0.5647)/(length(Pulseall)-1) 1-(Pulse-1)*(1-0)/(length(Pulseall)-1)], 'LineWidth', 1.5)
    end
    lgd=legend('1','2','3','4','5');
    lgd.Title.String = '# pulses';
    axis square
    title('Peak response all data')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');
    %for p2 plot the average
    p(1,2).select();
    hold on
    for Pulse=1:length(Pulseall)
        errorbar(AMPall,ratespiking_mean(:,Pulse),ratespiking_std(:,Pulse),'color',[0.0784+(Pulse-1)*(1-0.0784)/(length(Pulseall)-1) 0.5647+(Pulse-1)*(0.8431-0.5647)/(length(Pulseall)-1) 1-(Pulse-1)*(1-0)/(length(Pulseall)-1)], 'LineWidth', 1.5)
    end
         lgd=legend('1','2','3','4','5');
        lgd.Title.String = '# pulses';
        axis square
            title('All data Average 1-30ms after stim')
            ylabel('Sp/s')
            xlabel('Current \muA')
            set(gca,'TickDir','out');

figure;
    p = panel;
    p.pack(1,2)
    p(1,1).select();
    hold on
    for AMP=1:length(AMPall)
        errorbar(Pulseall,ratepeak(AMP,:),ratepeakstd(AMP,:),'Color',[0 AMP/length(AMPall) 1/AMP])
    end
    lgd=legend('2','5','6','8','10');
    lgd.Title.String = 'Current \muA';
    axis square
    title('All data Peak response')
    ylabel('Sp/s')
    xlabel('# Pulses')
    set(gca,'TickDir','out');
    %for p2 plot the average
    p(1,2).select();
    hold on
    for AMP=1:length(AMPall)
        errorbar(Pulseall,ratespiking_mean(AMP,:),ratespiking_std(AMP,:),'Color',[0 AMP/length(AMPall) 1/AMP])
    end

    lgd=legend('2','5','6','8','10');
    lgd.Title.String = 'Current \muA';
    axis square
            title('All data Average 1-30ms after stim')
            ylabel('Sp/s')
    xlabel('# Pulses')
            set(gca,'TickDir','out');

%% plot surpression now with all significant channels

ratespiking_suppression_mean=zeros(length(AMPall),length(Pulseall));
ratespiking_suppression_std=zeros(length(AMPall),length(Pulseall));
ratepeak_suppression=zeros(length(AMPall),length(Pulseall));
ratepeakstd_suppression=zeros(length(AMPall),length(Pulseall));
numsig_suppression=zeros(length(AMPall),5);
numelect_suppression=zeros(length(AMPall),5);
StimElectcount_suppression=zeros(length(AMPall),5);
%loop through amp and number of pulses then plot
ratespiking_suppression=ratespiking;
figure
    p = panel();
    p.pack(2, 3);
for AMP=1:length(AMPall)
    AMPcheck=['A' num2str(AMPall(AMP))];
    if AMP<4
        p(1,AMP).select();
    else
        p(2,AMP-3).select();
    end
    ax=p.de.axis;
    ax=ax(AMP);
    hold on
    
    for Pulse=1:length(Pulseall)
        Pulsecheck=['P' num2str(Pulseall(Pulse))];
        %make sure only responding elect included
       timesup=round(100+Pulse*3.333):175;
        if isfield(ratespiking_suppression.(AMPcheck),(Pulsecheck))
%         thresh=squeeze(mean(ratespiking_suppression.(AMPcheck).(Pulsecheck)(:,subselectbaseline),2)-(std(ratespiking_suppression.(AMPcheck).(Pulsecheck)(:,subselectbaseline),[],2).*SD_supression));
%         notsigsingle=squeeze(mean(ratespiking_suppression.(AMPcheck).(Pulsecheck)(:,timesup),2))>=thresh;
%         ratespiking_suppression.(AMPcheck).(Pulsecheck)(notsigsingle,:)=nan;
%         numsig_suppression(AMP,Pulse)=sum(~notsigsingle);
          
        numelect_suppression(AMP,Pulse)= sum(~isnan(ratespiking_suppression.(AMPcheck).(Pulsecheck)(:,1)));
        % Reshape the array into a 2D matrix with 128 elements per column
        reshapedArray = reshape(ratespiking_suppression.(AMPcheck).(Pulsecheck)(:,1), 64, []);
        % Count the number of non-NaN values in each column
        StimElectcount_suppression(AMP,Pulse) = sum(~isnan(sum(~isnan(reshapedArray))));
        maxdat=max(ratespiking_suppression.(AMPcheck).(Pulsecheck)(:,timesup),[],2);
        ratespiking_suppression_mean(AMP,Pulse)=nanmean(ratespiking_suppression.(AMPcheck).(Pulsecheck)(:,timesup),'all');
            %peakvals
        ratepeak_suppression(AMP,Pulse)=mean(maxdat,'omitnan');
        ratepeakstd_suppression(AMP,Pulse)=SEM(maxdat,0);
        stdshade(ratespiking_suppression.(AMPcheck).(Pulsecheck),0.2,[0.0784+(Pulse-1)*(1-0.0784)/(length(Pulseall)-1) 0.5647+(Pulse-1)*(0.8431-0.5647)/(length(Pulseall)-1) 1-(Pulse-1)*(1-0)/(length(Pulseall)-1)],[-90:300],1,ax)

        end
    end
    lgd=legend('1','2','3','4','5');
    lgd.Title.String = '# pulses';
    xlim([-90,90])
    title(['All data suppression' num2str(AMPall(AMP)) '\muA'])
    axis square
    ylabel('Sp/s')
    xlabel('Time (ms)')
    set(gca,'TickDir','out');
    ylim([-20 75])
    for j = 1:size(numelect_suppression,2)
        text(j*10-80, 65, num2str(numelect_suppression(AMP,j)), 'HorizontalAlignment', 'center');
    end

end

% p(2,3).select();
% hold on;
% 
% for Pulse=1:length(Pulseall)
% errorbar(AMPall,ratespiking_suppression_mean(:,Pulse),ratespiking_suppression_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
% end
%  lgd=legend('1','2','3','4','5');
% lgd.Title.String = '# pulses';
%     title('Average 1-40ms after stim')
%     ylabel('Sp/s')
%     xlabel('Current \muA')
%     set(gca,'TickDir','out');
    
    figure;
    p = panel;
    p.pack(1,2)
    p(1,1).select();
    hold on
    for Pulse=1:length(Pulseall)
        errorbar(AMPall,ratepeak_suppression(:,Pulse),ratepeakstd_suppression(:,Pulse),'color',[0.0784+(Pulse-1)*(1-0.0784)/(length(Pulseall)-1) 0.5647+(Pulse-1)*(0.8431-0.5647)/(length(Pulseall)-1) 1-(Pulse-1)*(1-0)/(length(Pulseall)-1)], 'LineWidth', 1.5)
    end
    lgd=legend('1','2','3','4','5');
    lgd.Title.String = '# pulses';
    axis square
    title('All data suppression Peak response')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');
    %for p2 plot the average
    p(1,2).select();
    hold on
    for Pulse=1:length(Pulseall)
        errorbar(AMPall,ratespiking_suppression_mean(:,Pulse),ratespiking_suppression_std(:,Pulse),'Color',[0.0784+(Pulse-1)*(1-0.0784)/(length(Pulseall)-1) 0.5647+(Pulse-1)*(0.8431-0.5647)/(length(Pulseall)-1) 1-(Pulse-1)*(1-0)/(length(Pulseall)-1)], 'LineWidth', 1.5)
    end
         lgd=legend('1','2','3','4','5');
        lgd.Title.String = '# pulses';
        axis square
            title('All data suppression Average 1-30ms after stim')
            ylabel('Sp/s')
            xlabel('Current \muA')
            set(gca,'TickDir','out');

figure;
    p = panel;
    p.pack(1,2)
    p(1,1).select();
    hold on
    for AMP=1:length(AMPall)
        errorbar(Pulseall,ratepeak_suppression(AMP,:),ratepeakstd_suppression(AMP,:),'Color',[0 AMP/length(AMPall) 1/AMP])
    end
    lgd=legend('2','5','6','8','10');
    lgd.Title.String = 'Current \muA';
    axis square
    title('All data suppression Peak response')
    ylabel('Sp/s')
    xlabel('# Pulses')
    set(gca,'TickDir','out');
    %for p2 plot the average
    p(1,2).select();
    hold on
    for AMP=1:length(AMPall)
        errorbar(Pulseall,ratespiking_suppression_mean(AMP,:),ratespiking_suppression_std(AMP,:),'Color',[0 AMP/length(AMPall) 1/AMP])
    end

    lgd=legend('2','5','6','8','10');
    lgd.Title.String = 'Current \muA';
    axis square
            title('All data suppression Average 1-30ms after stim')
            ylabel('Sp/s')
    xlabel('# Pulses')
            set(gca,'TickDir','out');



    %find nearest neighbour stim chn and add in the unique pulse data
    % use pensplit data to fill in missing data from the stim chn. FIrst ensure there are at least 3 different pulse types e.g. P1 P4 P5 or P2 P3 P4 then fill in the missing data
    % if there are only 2 pulse types then discard chanel from data
    % if there is only 1 pulse type then discard channel from data
    chnkeepmin_check=0;
   ratespiking_pensplit_stitch=[];
uniquepen=fields(ratespiking_pensplit);
for pen=1:length(uniquepen)
    uniqueamp=fields(ratespiking_pensplit.(uniquepen{pen}));
    for AMP=1:length(uniqueamp)
        uniquechn=fields(ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}));
        for chn=1:length(uniquechn)
            uniquepulse=fields(ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}));
            if length(uniquepulse)>2
                 ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn})=ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn});
                for Pulse=1:5
                    Pulsecheck=['P' num2str(Pulse)];
                    
                    %stitch all data from same stimchn first then deal with missing pulse
                    for pulsit=5:-1:Pulse+1
                        Pulsecheckstitch=['P' num2str(pulsit)];
                        if isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheck)) && isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheckstitch))
                        ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheckstitch)(:,94+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse))=ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheck)(:,94+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse));
                        else %find the closest match pulse
                            %get number from uniquechn instead of 'CHNxxx'.
                            chnnum_all_text=regexp(uniquechn,'\d+','match');
                            % Convert all chnnums to numerical value
                            chnnum_all  = cellfun(@str2double, chnnum_all_text);
                            chninterest=chnnum_all(chn);
                            chnnum_all_text=uniquechn;
                            chnnum_all_text(chn)=[];
                            chnnum_all(chn)=[];
                            % find nearest neighbour stim chn to chninterest
                            for i=1:size(chnnum_all)
                                [chnmin,idx] = min(abs(chninterest-chnnum_all));
                                if chnkeepmin_check<chnmin
                                chnkeepmin_check=chnmin;
                                end
                                if ~isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheck)) && isfield(ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(chnnum_all_text{idx}),(Pulsecheck)) 
                                    ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheck)=ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(chnnum_all_text{idx}).(Pulsecheck);
                                      break;
                                elseif isfield(ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(chnnum_all_text{idx}),(Pulsecheckstitch)) && ~isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheckstitch)) && isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheck))
                                    ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheckstitch)=ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(chnnum_all_text{idx}).(Pulsecheckstitch);
                                    ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheckstitch)(:,94+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse))=ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheck)(:,94+ceil(windowsize*(Pulse-1)):93+ceil(windowsize*Pulse));
                                    break;
                                elseif isfield(ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(chnnum_all_text{idx}),(Pulsecheckstitch)) && ~isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheckstitch)) && ~isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheck))
                                    ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheckstitch)=ratespiking_pensplit.(uniquepen{pen}).(uniqueamp{AMP}).(chnnum_all_text{idx}).(Pulsecheckstitch);
                                else
                                    chnnum_all_text(idx)=[];
                                    chnnum_all(idx)=[];
                                end
                            end
                            
                        end
                    end

                end 
            end
        end
    end
end
    %now combine and plot ratespiking_pensplit_stitch for each pulse and for each current
datallpenspit_stitch=[];
spatiallyorganisedarray=[];
    for pen=1:length(uniquepen)
        uniqueamp=fields(ratespiking_pensplit_stitch.(uniquepen{pen}));
        for AMP=1:length(uniqueamp)
            uniquechn=fields(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}));
            for chn=1:length(uniquechn)
                for Pulse=1:length(Pulseall)
                    Pulsecheck=['P' num2str(Pulseall(Pulse))];
                    if isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheck)) && isfield(datallpenspit_stitch,(AMPcheck)) && isfield(datallpenspit_stitch.(AMPcheck),(Pulsecheck))
                        datallpenspit_stitch.(uniqueamp{AMP}).(Pulsecheck)=cat(1,datallpenspit_stitch.(uniqueamp{AMP}).(Pulsecheck),ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheck));
                        spatiallyorganisedarray.(uniqueamp{AMP}).(Pulsecheck)=cat(3,spatiallyorganisedarray.(uniqueamp{AMP}).(Pulsecheck),plotheatmapdistdata(mean(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheck)(:,94:94+30),2), uniquechn{chn}));
                    elseif isfield(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}),(Pulsecheck))
                        datallpenspit_stitch.(uniqueamp{AMP}).(Pulsecheck)=ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheck);
                        spatiallyorganisedarray.(uniqueamp{AMP}).(Pulsecheck)=plotheatmapdistdata(mean(ratespiking_pensplit_stitch.(uniquepen{pen}).(uniqueamp{AMP}).(uniquechn{chn}).(Pulsecheck)(:,94:94+30),2), uniquechn{chn});
                    end
                    
                end
            end
        end
    end
    
    % now plot the stitched data
    figure
    p = panel();
    p.pack(2, 3);
    for AMP=1:length(AMPall)
        AMPcheck=['A' num2str(AMPall(AMP))];
        if AMP<4
            p(1,AMP).select();
        else
            p(2,AMP-3).select();
        end
        ax=p.de.axis;
        ax=ax(AMP);
        hold on
        for Pulse=1:length(Pulseall)
            Pulsecheck=['P' num2str(Pulseall(Pulse))];
            if isfield(datallpenspit_stitch,(AMPcheck)) && isfield(datallpenspit_stitch.(AMPcheck),(Pulsecheck))
                stdshade(datallpenspit_stitch.(AMPcheck).(Pulsecheck),0.2,[0.0784+(Pulse-1)*(1-0.0784)/(length(Pulseall)-1) 0.5647+(Pulse-1)*(0.8431-0.5647)/(length(Pulseall)-1) 1-(Pulse-1)*(1-0)/(length(Pulseall)-1)],[-90:300],1,ax)
            end
        end
        lgd=legend('1','2','3','4','5');
        lgd.Title.String = '# pulses';
        axis square
        xlim([-90,90])
        title(['All data stitch' num2str(AMPall(AMP)) '\muA'])
        ylabel('Sp/s')
        xlabel('Time (ms)')
        set(gca,'TickDir','out');
        ylim([0 20])
    end

%plot mean in third dimension of spatiallyorganisedarray
figure
    p = panel();
    p.pack(length(AMPall), length(Pulseall));
    for AMP=1:length(AMPall)
        AMPcheck=['A' num2str(AMPall(AMP))];
 
        for Pulse=1:length(Pulseall)
            p(AMP,Pulse).select();
            Pulsecheck=['P' num2str(Pulseall(Pulse))];
            if isfield(spatiallyorganisedarray,(AMPcheck)) && isfield(spatiallyorganisedarray.(AMPcheck),(Pulsecheck))
                imagesc(mean(spatiallyorganisedarray.(AMPcheck).(Pulsecheck),3,'omitnan'))
            end
            title([ num2str(AMPall(AMP)) '\muA' num2str(Pulseall(Pulse)) ' pulses'])
            axis square
            set(gca,'TickDir','out');
            
            clim([-5 5])
        end
        colorbar


    end


%     uniqueym=fields(Ratestitchspiking);
%     for itYM=1:length(uniqueym)
%         %iterate through Ratestitchspiking and average the data
%         uniqueamp=fields(Ratestitchspiking.(uniqueym{itYM}));
%         for itAMP=1:length(uniqueamp)
%             uniquechn=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}));
%             for itCHN=1:length(uniquechn)
%                 
%                 uniquepulse=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}));
%                 for itPulse=1:5
%                     Pulsecheck=['P' num2str(itPulse)];
%                     if ~isfield(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}),Pulsecheck) || size(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheck),2)~=391
%                         %get number from uniquechn instead of 'CHNxxx'.
%                         chnnum_all_text=regexp(uniquechn,'\d+','match');
%                         % Convert all chnnums to numerical value
%                         chnnum_all  = cellfun(@str2double, chnnum_all_text);
%                         chninterest=chnnum_all(itCHN);
%                         chnnum_all_text=uniquechn;
%                         chnnum_all_text(itCHN)=[];
%                         chnnum_all(itCHN)=[];
%                         % find nearest neighbour stim chn to chninterest
%                         for i=1:size(chnnum_all)
%                             [~,idx] = min(abs(chninterest-chnnum_all));
%                             if isfield(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}),Pulsecheck) && size(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck),2)==391
%                                 Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheck)=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck);
%                                 %add in pulse data from early pulses into late pulse
%                                 for pulsit=5:-1:itPulse
%                                     Pulsecheckstitch=['P' num2str(pulsit)];
%                                     if itPulse==pulsit && size(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheck),2)~=391
%                                         Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheckstitch)(:,94+ceil(windowsize*(itPulse-1)):391)=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck)(:,94+ceil(windowsize*(itPulse-1)):391);                                      
%                                     else
%                                         Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheckstitch)(:,94+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse))=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck)(:,94+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse));                                      
% %                                     else
% %                                         Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheckstitch)(:,91+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse))=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck)(chnrng,91+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse));
%                                     end
%                                 end
%                                 break;
%                             else
%                                 chnnum_all_text{idx}=[];
%                                 chnnum_all(idx)=nan;
%                                 
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end

    
    
    
%     
%     Ratestitchspikingmean=[];
%     for itYM=1:length(uniqueym)
%         %iterate through Ratestitchspiking and average the data
%         uniqueamp=fields(Ratestitchspiking.(uniqueym{itYM}));
%         for itAMP=1:length(uniqueamp)
%             uniquechn=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}));
%             for itCHN=1:length(uniquechn)
%                 uniquepulse=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}));
%                 for itPulse=1:length(uniquepulse)
%                     %put each unique combination of amp and pulse into a new array so that we can average them later ignoring uniqueym and chn
%                     if length(uniquepulse)>3 && size(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).P5,2)==391
%                         if ~isfield(Ratestitchspikingmean,(uniqueamp{itAMP})) || ~isfield(Ratestitchspikingmean.(uniqueamp{itAMP}),(uniquepulse{itPulse}))
%                             Ratestitchspikingmean.(uniqueamp{itAMP}).(uniquepulse{itPulse})=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(uniquepulse{itPulse});
%                         else
%                             Ratestitchspikingmean.(uniqueamp{itAMP}).(uniquepulse{itPulse})=cat(1,Ratestitchspikingmean.(uniqueamp{itAMP}).(uniquepulse{itPulse}),Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(uniquepulse{itPulse}));
%                         end
%                     end
% 
%                 end
%             end
%         end
%     
%     end
    %plot Ratestitchspikingmean for each current with the lines on the plot being the pulse data
%     figure;
%     p = panel();
%     p.pack(2, 3);
%     ratestitchmean=zeros(length(AMPall),length(Pulseall));
%     ratestitchstd=zeros(length(AMPall),length(Pulseall));
%     for AMP=1:length(AMPall)
%         AMPcheck=['A' num2str(AMPall(AMP))];
%         if AMP<4
%             p(1,AMP).select();
%         else
%             p(2,AMP-3).select();
%         end
%         ax=p.de.axis;
%         ax=ax(AMP);
%         hold on
%         for Pulse=1:length(Pulseall)
%             Pulsecheck=['P' num2str(Pulseall(Pulse))];
%             if isfield(Ratestitchspikingmean,(AMPcheck)) && isfield(Ratestitchspikingmean.(AMPcheck),(Pulsecheck))
%                 stdshade(Ratestitchspikingmean.(AMPcheck).(Pulsecheck),0.2,[0 Pulse/length(Pulseall) 1/Pulse],[-90:300],1,ax)
%                 ratestitchmean(AMP,Pulse)=nanmean(Ratestitchspikingmean.(AMPcheck).(Pulsecheck)(:,avgtime),'all');
%                 ratestitchstd(AMP,Pulse)=SEM(Ratestitchspikingmean.(AMPcheck).(Pulsecheck)(:,avgtime),0);
%             end
%         end
%         xlim([-90,90])
%         title([num2str(AMPall(AMP)) '\muA'])
%         ylabel('Sp/s')
%         xlabel('Time (ms)')
%         set(gca,'TickDir','out');
%     end
%     lgd=legend('1','2','3','4','5');
%     lgd.Title.String = '# pulses';
% p(2,3).select();
% hold on
% for Pulse=1:length(Pulseall)
%     errorbar(AMPall,ratestitchmean(:,Pulse),ratestitchstd(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
% 
% end
%  lgd=legend('1','2','3','4','5');
% lgd.Title.String = '# pulses';
%     title('Average 1-20ms after stim,Pulses added')
%     ylabel('Sp/s')
%     xlabel('Current \muA')
%     set(gca,'TickDir','out');

%% This is to plot the limited data of 2 monkeys and 2 penetrations with ALL pulse combinations
% plot the data for each pen if there are pulses 1-5 under the stim chn
% for each pen at each amp and each stimchn go in and check if pulse 1-5 are there
% if they are then add that data to a structure with amp and pulse as fields
% then plot the data
datincludestack=[];
stitchdata=[];
countit=0;
%initialise spikecountStitch for amplitudes
for AMP=1:length(AMPall)
    AMPcheck=['A' num2str(AMPall(AMP))];
        spikecountStitch.(AMPcheck)=nan(64,5,64*length(startpen));
        spikecountStitchsup.(AMPcheck)=nan(64,5,64*length(startpen));
end

for pen=1:length(startpen)
    pencheck=['Pen' num2str(pen)];
    if isfield(ratespiking_pensplit,(pencheck))
        for AMP=1:length(AMPall)
            AMPcheck=['A' num2str(AMPall(AMP))];
            if isfield(ratespiking_pensplit.(pencheck),(AMPcheck))
                for chn=65:128
                    stimchncheck=['SChn' num2str(chn)];
                    countit=(pen-1)*65+chn-64;
                    if isfield(ratespiking_pensplit.(pencheck).(AMPcheck),(stimchncheck))
                        %check if there are fields P2,P3,P4 and P5
                        if isfield(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck),'P1') && isfield(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck),'P2') && isfield(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck),'P3') && isfield(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck),'P4') && isfield(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck),'P5')
                            for Pulse=1:5
                                Pulsecheck=['P' num2str(Pulse)];
                                avgstitchdata.(AMPcheck).(Pulsecheck)=[];
                                %check significant channels in ratespiking_pensplit
                                % thresh=squeeze(mean(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)(:,subselectbaseline),2)+(std(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)(:,subselectbaseline),[],2).*SD_supression));
                                % notsigsingle=squeeze(mean(ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)(:,timepeak),2))<=thresh;    
                                % ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)(notsigsingle,:)=nan;
                                if isfield(datincludestack,(AMPcheck)) && isfield(datincludestack.(AMPcheck),(Pulsecheck))
                                    datincludestack.(AMPcheck).(Pulsecheck)=cat(1,datincludestack.(AMPcheck).(Pulsecheck),ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck));
                                else
                                    datincludestack.(AMPcheck).(Pulsecheck)=ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck);
                                end
                                for pulseit=5:-1:Pulse
                                     Pulsecheck2=['P' num2str(pulseit)];
                                     if pulseit==Pulse
                                        stitchdata.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck2)(:,94+ceil(windowsize*(Pulse-1)):391)=ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)(:,94+ceil(windowsize*(Pulse-1)):391);
                                        stitchdata.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck2)(:,1:90)=ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)(:,1:90);
                                     elseif pulseit>Pulse
                                         stitchdata.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck2)(:, 94+ceil(windowsize*(Pulse-1)):94+ceil(windowsize*(Pulse)))=ratespiking_pensplit.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)(:, 94+ceil(windowsize*(Pulse-1)):94+ceil(windowsize*(Pulse)));
                                     end
                                end
                            end
                            for pulse=1:5
                                pcheck=(['P' num2str(pulse)]);
                                stitchdata.(pencheck).(AMPcheck).(stimchncheck).(pcheck)(sum(stitchdata.(pencheck).(AMPcheck).(stimchncheck).(pcheck)(:,94:181),2)==0,:)=nan;
                            %add the numbers between columns 94 and 120 of the stitchdata then average the counts for each row, stim channel and pulse
                            
                            spikecountStitch.(AMPcheck)(1:64,pulse,countit)=sum(stitchdata.(pencheck).(AMPcheck).(stimchncheck).(pcheck)(:,94:120),2);%94:120
                            spikecountStitchsup.(AMPcheck)(1:64,pulse,countit)=sum(stitchdata.(pencheck).(AMPcheck).(stimchncheck).(pcheck)(:,121:185),2);

                            end
                        end
                    end
                end
            end
        end
    end
end
% plot spikecountstitch average for each pulse on a figure for each amplitude
figure;
p=panel();
p.pack(2,2);
spikcount=nan(5,5);
spikecountsem=nan(5,5);
supspikecount=nan(5,5);
spikecountsemsup=nan(5,5);
dataforanova=[];
dataforanovasup=[];
F1=[];
F2=[];
numelect=nan(5,5);
for AMP=1:length(AMPall)
    AMPcheck=['A' num2str(AMPall(AMP))];

        if isfield(spikecountStitch,(AMPcheck)) 
            p(2,1).select();
            hold on;
            dat=reshape(permute(spikecountStitch.(AMPcheck),[2,1,3]),[5,64*size(spikecountStitch.(AMPcheck),3)]);
            dataforanova=[dataforanova, dat];
            F1=[F1,ones(size(dat)).*AMPall(AMP)];
            F2=[F2,repmat([1; 2; 3; 4; 5],[1,size(dat,2)])];
            spikcount(AMP,:)=mean(dat,2,'omitnan');
            spikecountsem(AMP,:)=SEM(dat,1);
           % countelect(AMP,:)=sum(~isnan(dat),2);
            errorbar(1:5,spikcount(AMP,:),spikecountsem(AMP,:),'Color',map3(AMP,:),'LineWidth', 1.5) 
            p(2,2).select();
            hold on;
            dat2=reshape(permute(spikecountStitchsup.(AMPcheck),[2,1,3]),[5,64*size(spikecountStitch.(AMPcheck),3)]);
            errorbar(1:5,mean(dat2,2,'omitnan'),SEM(dat2,1),'Color',map3(AMP,:),'LineWidth', 1.5) 
             dataforanovasup=[dataforanovasup, dat2];
             supspikecount(AMP,:)=mean(dat2,2,'omitnan');
              spikecountsemsup(AMP,:)=SEM(dat2,1);
              numelect(AMP,:)=sum(~isnan(dat),2);
        end

end

p(2,1).select();
    title('Select data Stitched spike count response 4-30ms after stim')
    ylabel('Sp/s')
    xlabel('# Pulses')
    set(gca,'TickDir','out');
    axis square
p(2,2).select();
    title('Select data Stitched spike count response 31-85ms after stim')
    ylabel('Sp/s')
    xlabel('# Pulses')
    set(gca,'TickDir','out');
    lgd=legend('2','5','6','8','10');
    lgd.Title.String = 'Current \muA';
    axis square
    
for pulse=1:5
     p(1,1).select();
            hold on;
 errorbar(AMPall,spikcount(:,pulse),spikecountsem(:,pulse),'Color',map4(pulse,:), 'LineWidth', 1.5)
 p(1,2).select();
 hold on;
  errorbar(AMPall,supspikecount(:,pulse),spikecountsemsup(:,pulse),'Color',map4(pulse,:),'LineWidth', 1.5)
end

p(1,1).select();
    title('Select data Stitched spike count response 4-30ms after stim')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');
    axis square
p(1,2).select();
    title('Select data Stitched spike count response 31-85ms after stim')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');
    lgd=legend('1','2','3','4','5');
    lgd.Title.String = '# Pulses';
    axis square



[p,tbl,stats]=anovan(dataforanova(:),{F1(:),F2(:)},"Varnames",["Current","Pulse"]);
axis square
[results,tbl,h,gnames] = multcompare(stats,"Dimension",[1 2]);
title('4-30ms stitched data spike count')
set(gca,'TickDir','out');
fontname('Times New Roman')
axis square
x=categorical(gnames);
x=reordercats(x,gnames);
nm=char(x);
lblnm=[nm(:,9) repmat(['\muA'],size(nm,1),1)];

     bar_werror(lblnm,spikcount(:),spikecountsem(:))
     axis square
     fontname('Times New Roman')
   

  grid_significance(x,results)
axis square
%count how many significantly different comparisons there are with changing current. gnames has each condition listed as 'Current#Pulse#' and results has each comparison between conditions
%listed as [condition1 condition2 lowerCI upperCI pvalue]
Currentsamediff=[0 0;0 0];
for i=1:length(results)
    if gnames{results(i,1)}(9)==gnames{results(i,2)}(9)
        Currentsamediff(2,1)=Currentsamediff(2,1)+1;
        if results(i,6)<0.05
            Currentsamediff(1,1)=Currentsamediff(1,1)+1;
        end
    end
end
%do the same for different currents but same pulse count
for i=1:length(results)
    if (gnames{results(i,1)}(9)~=gnames{results(i,2)}(9)) && (gnames{results(i,1)}(end)==gnames{results(i,2)}(end))
        Currentsamediff(2,2)=Currentsamediff(2,2)+1;
        if results(i,6)<0.05
            Currentsamediff(1,2)=Currentsamediff(1,2)+1;
        end
    end
end
%now do the same for pulses
Pulsesamediff=[0 0;0 0];
for i=1:length(results)
    if gnames{results(i,1)}(end)==gnames{results(i,2)}(end)
        Pulsesamediff(2,1)=Pulsesamediff(2,1)+1;
        if results(i,6)<0.05
            Pulsesamediff(1,1)=Pulsesamediff(1,1)+1;
        end
    end
end
%do the same for different pulses
for i=1:length(results)
    if gnames{results(i,1)}(end)~=gnames{results(i,2)}(end)
        Pulsesamediff(2,2)=Pulsesamediff(2,2)+1;
        if results(i,6)<0.05
            Pulsesamediff(1,2)=Pulsesamediff(1,2)+1;
        end
    end
end

figure;
[p,tbl,stats]=anovan(dataforanovasup(:),{F1(:),F2(:)},"Varnames",["Current","Pulse"]);
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);
title('30-75ms stitched data spike count')
x=categorical(gnames);
x=reordercats(x,gnames);
nm=char(x);
lblnm=[nm(:,9) repmat(['\muA'],size(nm,1),1)];
        bar_werror(lblnm,supspikecount(:),spikecountsemsup(:))
    ylabel('Spike count 30-75ms after stimulus onset');
    axis square
    grid_significance(x,results)
    axis square
%count how many significantly different comparisons there are with changing current. gnames has each condition listed as 'Current#Pulse#' and results has each comparison between conditions
%listed as [condition1 condition2 lowerCI upperCI pvalue]
Currentsamediff=[0 0;0 0];
for i=1:length(results)
    if gnames{results(i,1)}(9)==gnames{results(i,2)}(9)
        Currentsamediff(2,1)=Currentsamediff(2,1)+1;
        if results(i,6)<0.05
            Currentsamediff(1,1)=Currentsamediff(1,1)+1;
        end
    end
end
%do the same for different currents
for i=1:length(results)
    if gnames{results(i,1)}(9)~=gnames{results(i,2)}(9) && (gnames{results(i,1)}(end)==gnames{results(i,2)}(end))
        Currentsamediff(2,2)=Currentsamediff(2,2)+1;
        if results(i,6)<0.05
            Currentsamediff(1,2)=Currentsamediff(1,2)+1;
        end
    end
end
%now do the same for pulses
Pulsesamediff=[0 0;0 0];
for i=1:length(results)
    if gnames{results(i,1)}(end)==gnames{results(i,2)}(end)
        Pulsesamediff(2,1)=Pulsesamediff(2,1)+1;
        if results(i,6)<0.05
            Pulsesamediff(1,1)=Pulsesamediff(1,1)+1;
        end
    end
end
%do the same for different pulses
for i=1:length(results)
    if gnames{results(i,1)}(end)~=gnames{results(i,2)}(end) && gnames{results(i,1)}(9)==gnames{results(i,2)}(9)
        Pulsesamediff(2,2)=Pulsesamediff(2,2)+1;
        if results(i,6)<0.05
            Pulsesamediff(1,2)=Pulsesamediff(1,2)+1;
        end
    end
end

    %collapsed data
    pulcol=nan(5,3);
     curcol=nan(5,3);
    for amp=1:length(AMPall)
    curcol(amp,1)=mean(dataforanova((F1==AMPall(amp) & F2==Pulseall(1))),'omitnan');
    curcol(amp,2)=SEM(dataforanova((F1==AMPall(amp) & F2==Pulseall(1))),0);
     curcol(amp,3)=sum(~isnan(dataforanova((F1==AMPall(amp) & F2==Pulseall(1)))));
    pulcol(amp,1)=mean(dataforanova(F1==AMPall(5) & F2==Pulseall(amp)),'omitnan');
    pulcol(amp,2)=SEM(dataforanova(F1==AMPall(5) & F2==Pulseall(amp)),0);
     pulcol(amp,3)=sum(~isnan(dataforanova(F1==AMPall(5) & F2==Pulseall(amp))));
    end
    figure
     errorbar(AMPall,curcol(:,1),curcol(:,2),'k')
     xlabel('Current (\muA)')
     ylabel('Spike count 4-10ms')
      set(gca,'TickDir','out');
      text(3, 60, [num2str(min(curcol(:,3))) '-' num2str(max(curcol(:,3)))]) 
      xlim([1 10])
        axis square
     figure
     errorbar(Pulseall,pulcol(:,1),pulcol(:,2),'k')
      xlabel('# Pulses')
       set(gca,'TickDir','out');
     ylabel('Spike count 4-30ms')
    text(1.5, 25, [num2str(min(pulcol(:,3))) '-' num2str(max(pulcol(:,3)))]) 
    xlim([0.5 5])
    axis square

    

    chargeinjection=AMPall.*Pulseall';
    figure;
    %make scatter plot with errorbars using SEM
    x=chargeinjection(:);
    [x, sortIdx] = sort(x, 'ascend');
    y=spikcount(:);
    y=y(sortIdx);
  
    errorbar(x,y,spikecountsem(:),'k*')
%       [uniqueVals, ~, idx] = unique(x);
%       meanVals=zeros(length(uniqueVals),1);
%       for i = 1:length(uniqueVals)
%           meanVals(i) = mean(y(idx == i));
%       end
%     y=meanVals;
%     x=uniqueVals;


% Define the Naka-Rushton function
nakaRushton = fittype('Rmax * (x^n) / (C^n + x^n)', ...
                      'independent', 'x', ...
                      'coefficients', {'Rmax', 'C', 'n'});

% Set initial parameter estimates
initialGuess = [max(y), mean(x), 2];

% Fit the model
[fitresult, gof] = fit(x, y, nakaRushton, 'StartPoint', initialGuess);

% Generate fitted values
x_fit = linspace(min(x), max(x), 100);
y_fit = feval(fitresult, x_fit);
% Plot the Naka-Rushton fit
hold on
plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

%     p = polyfit(x, y,3);
%     y_fit = polyval(p, x);
%     hold on;
%     plot(x, y_fit, 'r-', 'LineWidth', 2);
    title('Select data Stitched data spike count (4-30ms)')
    axis square
    xlabel('Total charge injection (\muA)')
    ylabel('Spike count')
     set(gca,'TickDir','out');



    %plot sititchdata time vs spiking for each current and pulse
    for AMP=1:length(AMPall)
        AMPcheck=['A' num2str(AMPall(AMP))];
        for Pulse=1:5
            Pulsecheck=['P' num2str(Pulse)];
            for pen=1:6
                pencheck=['Pen' num2str(pen)];
                for stimchn=65:128
                    stimchncheck=['SChn' num2str(stimchn)];
                    if isfield(stitchdata,(pencheck)) && isfield(stitchdata.(pencheck).(AMPcheck),(stimchncheck))
                        
                        avgstitchdata.(AMPcheck).(Pulsecheck)=[avgstitchdata.(AMPcheck).(Pulsecheck); stitchdata.(pencheck).(AMPcheck).(stimchncheck).(Pulsecheck)];
                        
                    end
                end
            end
        end
    end
%plot the average of the stitched data
figure
p = panel();
p.pack(2, 3);
significantstitch=nan(length(Pulseall),length(AMPall));
datascatterplot=nan(length(Pulseall),length(AMPall),size(avgstitchdata.A2.P1,1));
datatimerecove=nan(length(Pulseall),length(AMPall),size(avgstitchdata.A2.P1,1));
datasupp=nan(length(Pulseall),length(AMPall),size(avgstitchdata.A2.P1,1));
for AMP=1:length(AMPall)
    AMPcheck=['A' num2str(AMPall(AMP))];
    if AMP<4
        p(1,AMP).select();
    else
        p(2,AMP-3).select();
    end
    ax=p.de.axis;
    ax=ax(AMP);
     yline(0,'color','k')
    hold on
    for Pulse=1:5
        Pulsecheck=['P' num2str(Pulse)];
            stdshade(avgstitchdata.(AMPcheck).(Pulsecheck),0.2,map4(Pulse,:),[-90:300],5,ax)
            plot([3.333*(Pulse-1) 3.333*(Pulse-1)],[-3 -2.5],'color',map4(Pulse,:), 'LineWidth', 1.5)
            significantstitch(AMP,Pulse)=sum(~isnan(avgstitchdata.(AMPcheck).(Pulsecheck)(:,1)));
            datascatterplot(AMP,Pulse,:)=mean(avgstitchdata.(AMPcheck).(Pulsecheck)(:,94:180),2);
            dat=movmean(avgstitchdata.(AMPcheck).(Pulsecheck),5,2);
            for i=1:size(avgstitchdata.(AMPcheck).(Pulsecheck),2)
                datasupp(AMP,Pulse,i)=mean(avgstitchdata.(AMPcheck).(Pulsecheck)(i,120:180),2);
                
                if ~isnan(mean(dat(i,1:90),2))
                    %need to check when the data returns to baseline and remains there for 10samples
                    %find where data is within 1SD of baseline for 10 samples
                    baselinethresh(1)=mean(dat(i,1:90),2)-std(dat(i,1:90),[],2)./2;
                    baselinethresh(2)=mean(dat(i,1:90),2)+std(dat(i,1:90),[],2)./2;
                    withinbounds=dat(i,ceil(windowsize*(Pulse-1))+90:180)>=baselinethresh(1) & dat(i,ceil(windowsize*(Pulse-1))+90:180)<=baselinethresh(2);
                    %find where the consecutive samples are greater than 9 or equal to 10
                    consecutive_sum = conv(double(withinbounds), ones(1, 10), 'valid');
                    if any(consecutive_sum==10)
                        datatimerecove(AMP,Pulse,i)=find(consecutive_sum==10,1,'first');
                    else
                        datatimerecove(AMP,Pulse,i)=90;
                    end
                    %datatimerecove(AMP,Pulse,i)=find(dat(i,91:180)>=mean(dat(i,1:89),2),1,'first')+29;

                    %datatimerecove(AMP,Pulse,i)=find(avgstitchdata.(AMPcheck).(Pulsecheck)(i,120:180)>=mean(avgstitchdata.(AMPcheck).(Pulsecheck)(i,1:90),2),1,'first')+29;
                else
                    datatimerecove(AMP,Pulse,i)=nan;
                end
            end
    end
    xlim([0,90])
    title(['Select data ' num2str(AMPall(AMP)) '\muA'])
    axis square
    ylabel('Firing rate (Sp/s)')
    xlabel('Time (ms)')
    set(gca,'TickDir','out');
   
    ylim([-3 10])
end
lgd=legend('1','2','3','4','5');
lgd.Title.String = '# pulses';
% plot averaged stitched data using a heatmap with total charge injection on the y axis (current x pulses) and time on x axis from -90:90. Z axis is spiking rate
figure
iter=0;
chargeinjection=nan(length(AMPall)*length(Pulseall)*size(avgstitchdata.A2.P1,1),1);
datacharginjection=nan(391,length(AMPall)*length(Pulseall)*size(avgstitchdata.A2.P1,1));
  stdchargeinjection=datacharginjection;
for AMP=1:length(AMPall)
    AMPcheck=['A' num2str(AMPall(AMP))];
    for Pulse=1:5
        Pulsecheck=['P' num2str(Pulse)];
        iter=iter+1;
       chargeinjection(iter)=AMPall(AMP)*Pulse;
       datacharginjection(:,iter)=mean(avgstitchdata.(AMPcheck).(Pulsecheck),1,'omitnan');
       stdchargeinjection(:,iter)=SEM(avgstitchdata.(AMPcheck).(Pulsecheck),0);
       

    end
end
%then average all channels with same charge injection
chargeinjectionu=unique(chargeinjection);
chargeinjectionu(isnan(chargeinjectionu))=[];
datacharginjectionmean=nan(391,length(chargeinjectionu)); 
datasamechargeinjection=nan(391,5,length(AMPall)*length(Pulseall));   
stdchargeinjectionsame=datasamechargeinjection;
for i=1:length(chargeinjectionu)
    if ~isnan(chargeinjectionu(i))
        datacharginjectionmean(:,i)=mean(datacharginjection(:,chargeinjection==chargeinjectionu(i)),2);
        %pull those with the same charge injection but different pulse and current values and plot them
        datasamechargeinjection(:,1:sum(chargeinjection==chargeinjectionu(i)),i)=datacharginjection(:,chargeinjection==chargeinjectionu(i));
        stdchargeinjectionsame(:,1:sum(chargeinjection==chargeinjectionu(i)),i)=stdchargeinjection(:,chargeinjection==chargeinjectionu(i));
    end
end
%plot the data
figure
chargey=repmat(chargeinjectionu',[181,1]);
timex=repmat([-90:90]',[1,17]);
surf(timex,chargey,datacharginjectionmean(1:181,:))    
shading flat
ylabel('Current (\muA)')
xlabel('Time (ms)')
zlabel('Sp/s')
title('Select data Stitched data')
axis square
%make axis tight and adjust to display in x and y with a colourbar
view(2)
c=colorbar;
c.Label.String='Sp/s';

%plot the data for each charge injection on a heatmap only if there is more than one column of data
chargeinjectplot=[6 10 30 40];

for i=1:length(chargeinjectionu)
    if sum(~isnan(datasamechargeinjection(1,:,i)))>1
        figure
        surf(timex(:,1:5)',datasamechargeinjection(1:181,:,i)')
        title(['Charge injection ' num2str(chargeinjectionu(i))])
        axis square
        xlabel('Time (ms)')
        zlabel('Sp/s')
        axis tight
        view(2)
        shading flat
        colorbar
        %y labels will be charge injection divided by AMPall and only those that are whole numbers <5 e.g. 40 is 5*8 and 4*10 so 5p8ua and 4p10ua will be labels
        
        possiblelabels=chargeinjectionu(i)./AMPall';
        chargamp=AMPall(possiblelabels<=5);
        possiblelabels=possiblelabels(possiblelabels<=5);
        chargamp=chargamp(possiblelabels==round(possiblelabels));
        possiblelabels=possiblelabels(possiblelabels==round(possiblelabels));
        yticks(1.5:1:length(possiblelabels)+0.5)
        yticklabels([strcat(num2str(possiblelabels),{'p'},num2str(chargamp'),{'\muA'})])
        ylim([1 length(possiblelabels)+1])

        if any(chargeinjectionu(i)==chargeinjectplot)
            figure(25)
            ax=subplot(2,2,find(chargeinjectionu(i)==chargeinjectplot));
            hold on;
            
            for j=1:length(possiblelabels)
            %errorbar(timex(:,1),datasamechargeinjection(1:181,j,i),stdchargeinjectionsame(1:181,j,i))
            stdshade_errorcalc(datasamechargeinjection(1:181,j,i)',stdchargeinjectionsame(1:181,j,i)',0.2,map4(possiblelabels(j),:),timex(:,1),5,ax);
            title(['Charge injection ' num2str(chargeinjectionu(i))])
            axis square
            xlabel('Time (ms)')
            ylabel('Firing rate (Sp/s)')
            set(gca,'TickDir','out');
            ylim([-3 10])
            xlim([0 90])
            
            end
            legend([strcat(num2str(possiblelabels),{'p'},num2str(chargamp'),{'\muA'})])
            yline(0,'color','k')
            %xline(0,'r')
        end
    end
end






%recovery time from 30ms for select data
figure
xamp=repmat(AMPall',[1,5]);
ypulse=repmat(Pulseall,[5,1]);
zdat=mean(datatimerecove,3,'omitnan');
zsem=std(datatimerecove,0,3,'omitnan')/sqrt(size(datatimerecove,3));
surf(xamp,ypulse,zdat)%% this isn't working
%pcolor(AMPall,Pulseall,zdat)
xlabel('Current (\muA)')
ylabel('# Pulses')
c=colorbar;
c.Label.String='Time (ms)';
title('Time to return to baseline following suppresion')
axis square

figure
p=panel();
p.pack(1,2);
p(1,1).select();
colormap(map4)
title('Time to return to baseline following suppresion')
hold on
[f,~]=fit(xamp(:),zdat(:),'poly1');
for i=1:5
errorbar(xamp(:,i),zdat(:,i),zsem(:,i),'Color',map4(i,:),"LineStyle",'none',"Marker",'*');
end
plot(f,'k')
xlabel('Current (\muA)')
ylabel('Time (ms)')
set(gca,'TickDir','out');
axis square
xlim([1 10])
p(1,2).select();
title('Time to return to baseline following suppresion')
hold on
colormap(map3)
[f,~]=fit(ypulse(:),zdat(:),'poly1');
for i=1:5
errorbar(ypulse(i,:),zdat(i,:),zsem(i,:),'Color',map3(i,:),"LineStyle",'none',"Marker",'*');
end
plot(f,'k')
xlabel('# Pulses')
ylabel('Time (ms)')
set(gca,'TickDir','out');
axis square
xlim([0 5])
mdl=fitlm([ypulse(:),xamp(:)],zdat(:));%% for significance testing

%avg suppression
% figure
% xamp=repmat(AMPall',[1,5]);
% ypulse=repmat(Pulseall,[5,1]);
% zdat=mean(datasupp,3,'omitnan');
% scatter3(xamp(:),ypulse(:),zdat(:),100,zdat(:),'filled')
% xlabel('Current (\muA)')
% ylabel('# Pulses')
% c=colorbar;
% c.Label.String='Sp/s';
% % c.Limits=[min(zdat(:)) max(zdat(:))];
% % caxis([min(zdat) max(zdat)])
% title('Mean supp sp/s')



%now plot datincludestack also make sure the data is significant
datincludestack_supression=datincludestack;
ratespiking_mean=zeros(length(AMPall),length(Pulseall));    
ratespiking_std=zeros(length(AMPall),length(Pulseall));
ratepeak=zeros(length(AMPall),length(Pulseall));
ratepeakstd=zeros(length(AMPall),length(Pulseall));
numelect=zeros(length(AMPall)+1,5+1);
numelect(1,2:end)=Pulseall;
numelect(2:end,1)=AMPall';
figure;
p = panel();
p.pack(2, 3);
for AMP=1:length(AMPall)
    AMPcheck=['A' num2str(AMPall(AMP))];
    if AMP<4
        p(1,AMP).select();
    else
        p(2,AMP-3).select();
    end
    ax=p.de.axis;
    ax=ax(AMP);
    hold on
    for Pulse=1:5
        Pulsecheck=['P' num2str(Pulse)];
        timepeak=round(90+Pulse*3.333):round(120+Pulse*3.333);
        if isfield(datincludestack,(AMPcheck)) && isfield(datincludestack.(AMPcheck),(Pulsecheck))
%             thresh = squeeze(mean(datincludestack.(AMPcheck).(Pulsecheck)(:,subselectbaseline),2)+(std(datincludestack.(AMPcheck).(Pulsecheck)(:,subselectbaseline),[],2).*3));
%             notsigsingle = squeeze(mean(datincludestack.(AMPcheck).(Pulsecheck)(:,timepeak),2))<=thresh;
%             datincludestack.(AMPcheck).(Pulsecheck)(notsigsingle,:)=nan;
            numelect(AMP+1,Pulse+1)= sum(~isnan(datincludestack.(AMPcheck).(Pulsecheck)(:,1)));
            length(notsigsingle)-sum(notsigsingle)
            ratespiking_mean(AMP,Pulse)=nanmean(datincludestack.(AMPcheck).(Pulsecheck)(:,timepeak),'all');%nanmean(maxdat(maxdat~=0),'all');%nanmean(ratespiking.(AMPcheck).(Pulsecheck)(:,avgtime),'all');
            ratespiking_std(AMP,Pulse)=SEM(datincludestack.(AMPcheck).(Pulsecheck)(:,timepeak),0);
            %peakvals
            maxdat=max(datincludestack.(AMPcheck).(Pulsecheck)(:,round(90+Pulse*3.333):round(95+Pulse*3.333)),[],2);
            ratepeak(AMP,Pulse)=mean(maxdat,'omitnan');
            ratepeakstd(AMP,Pulse)=SEM(maxdat,0);
            stdshade(datincludestack.(AMPcheck).(Pulsecheck),0.2,map4(Pulse,:),[-90:300], 1,ax);
        end
    end
    xlim([-90,90])
    title(['Select data NOT stitch ' num2str(AMPall(AMP)) '\muA'])
    axis square
    ylabel('Sp/s')
    xlabel('Time (ms)')
    set(gca,'TickDir','out');
    ylim([0 95])
    for j = 2:size(numelect,2)
        text(j*10-80, 65, num2str(numelect(AMP+1,j)), 'HorizontalAlignment', 'center');
    end

end
% p(2,3).select();
% hold on
% for Pulse=1:5
%     errorbar(AMPall,ratespiking_mean(:,Pulse),ratespiking_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
% end
% lgd=legend('2','3','4','5');
% lgd.Title.String = '# pulses';
% title('Average 1-40ms after stim')
% ylabel('Sp/s')
% xlabel('Current \muA')
% set(gca,'TickDir','out');

%plot the peak data and average
figure;
p = panel;
p.pack(1,2)
p(1,1).select();
hold on
for Pulse=1:5
    errorbar(AMPall,ratepeak(:,Pulse),ratepeakstd(:,Pulse),'color',map4(Pulse,:), 'LineWidth', 1.5)
end
lgd=legend('2','3','4','5');
lgd.Title.String = '# pulses';
title('Select data NOT stitch Peak response')
ylabel('Sp/s')
xlabel('Current \muA')
set(gca,'TickDir','out');
axis square
%for p2 plot the average
p(1,2).select();
hold on
for Pulse=1:5
    errorbar(AMPall,ratespiking_mean(:,Pulse),ratespiking_std(:,Pulse),'color',map4(Pulse,:), 'LineWidth', 1.5)
end

lgd=legend('2','3','4','5');
lgd.Title.String = '# pulses';
title('Select data NOT stitch Average 1-30ms after stim')
ylabel('Sp/s')
xlabel('Current \muA')
set(gca,'TickDir','out');
axis square

%% plot against Pulse on x axis and current as different lines

figure;
p = panel;
p.pack(1,2)
p(1,1).select();
hold on
for AMP=1:length(AMPall)
    errorbar(Pulseall,ratepeak(AMP,:),ratepeakstd(AMP,:),'Color',map3(AMP,:))
end
lgd=legend('2','5','6','8','10');
lgd.Title.String = 'Current \muA';
title('Select data NOT stitch Peak response')
axis square
ylabel('Sp/s')
xlabel('# Pulses')
set(gca,'TickDir','out');
%for p2 plot the average
p(1,2).select();
hold on
for AMP=1:length(AMPall)
    errorbar(Pulseall,ratespiking_mean(AMP,:),ratespiking_std(AMP,:),'Color',map3(AMP,:))
end
lgd=legend('2','5','6','8','10');
lgd.Title.String = 'Current \muA';
axis square
title('Select data NOT stitch Average 1-30ms after stim')
ylabel('Sp/s')
xlabel('# Pulses')
set(gca,'TickDir','out');
end
%%
% %suppression data
% ratespiking_suppression_mean=zeros(length(AMPall),length(Pulseall));
% ratespiking_suppression_std=zeros(length(AMPall),length(Pulseall));
% ratepeak_suppression=zeros(length(AMPall),length(Pulseall));
% ratepeakstd_suppression=zeros(length(AMPall),length(Pulseall));
% numelect=zeros(length(AMPall)+1,5+1);
% numelect(1,2:end)=Pulseall;
% numelect(2:end,1)=AMPall';
% 
% figure;
% p = panel();
% p.pack(2, 3);
% for AMP=1:length(AMPall)
%     AMPcheck=['A' num2str(AMPall(AMP))];
%     if AMP<4
%         p(1,AMP).select();
%     else
%         p(2,AMP-3).select();
%     end
%     ax=p.de.axis;
%     ax=ax(AMP);
%     hold on
%     for Pulse=1:5
%         Pulsecheck=['P' num2str(Pulse)];
%         timesup=round(100+Pulse*3.333):175;
%         if isfield(datincludestack_supression,(AMPcheck)) && isfield(datincludestack_supression.(AMPcheck),(Pulsecheck))
%             thresh = squeeze(mean(datincludestack_supression.(AMPcheck).(Pulsecheck)(:,subselectbaseline),2)-(std(datincludestack_supression.(AMPcheck).(Pulsecheck)(:,subselectbaseline),[],2).*SD_supression));
%             notsigsingle = squeeze(mean(datincludestack_supression.(AMPcheck).(Pulsecheck)(:,timesup),2))>=thresh;
%             datincludestack_supression.(AMPcheck).(Pulsecheck)(notsigsingle,:)=nan;
%             numelect(AMP+1,Pulse+1)= sum(~isnan(datincludestack_supression.(AMPcheck).(Pulsecheck)(:,1)));
%             length(notsigsingle)-sum(notsigsingle)
%             ratespiking_suppression_mean(AMP,Pulse)=nanmean(datincludestack_supression.(AMPcheck).(Pulsecheck)(:,timesup),'all');
%             ratespiking_suppression_std(AMP,Pulse)=SEM(datincludestack_supression.(AMPcheck).(Pulsecheck)(:,timesup),0);
%             %peakvals
%             maxdat=max(datincludestack_supression.(AMPcheck).(Pulsecheck)(:,timesup),[],2);
%             ratepeak_suppression(AMP,Pulse)=mean(maxdat,'omitnan');
%             ratepeakstd_suppression(AMP,Pulse)=SEM(maxdat,0);
%             stdshade(datincludestack_supression.(AMPcheck).(Pulsecheck),0.2,[0 Pulse/length(Pulseall) 1/Pulse],[-90:300],1,ax);
%         end
%     end
%     xlim([-90,90])
%     title(['Select data NOT stitch suppression' num2str(AMPall(AMP)) '\muA'])
%     ylabel('Sp/s')
%     xlabel('Time (ms)')
%     set(gca,'TickDir','out');
%     ylim([-20 75])
%     for j = 2:size(numelect,2)
%         text(j*10-80, 65, num2str(numelect(AMP+1,j)), 'HorizontalAlignment', 'center');
%     end
% end
% % p(2,3).select();
% % hold on
% % for Pulse=1:5
% %     errorbar(AMPall,ratespiking_suppression_mean(:,Pulse),ratespiking_suppression_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
% % end
% % lgd=legend('2','3','4','5');
% % lgd.Title.String = '# pulses';
% % title('Average 1-40ms after stim')
% % ylabel('Sp/s')
% % xlabel('Current \muA')
% % set(gca,'TickDir','out');
% 
% %plot the peak data and average
% figure;
% p = panel;
% p.pack(1,2)
% p(1,1).select();
% hold on
% for Pulse=1:5
%     errorbar(AMPall,ratepeak_suppression(:,Pulse),ratepeakstd_suppression(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
% end
% lgd=legend('2','3','4','5');
% lgd.Title.String = '# pulses';
% title('Select data NOT stitch suppression Peak response')
% ylabel('Sp/s')
% xlabel('Current \muA')
% set(gca,'TickDir','out');
% %for p2 plot the average
% p(1,2).select();
% hold on
% for Pulse=1:5
%     errorbar(AMPall,ratespiking_suppression_mean(:,Pulse),ratespiking_suppression_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
% end
% lgd=legend('2','3','4','5');
% lgd.Title.String = '# pulses';
% title('Select data NOT stitch suppression Average 1-30ms after stim')
% ylabel('Sp/s')
% xlabel('Current \muA')
% set(gca,'TickDir','out');
% 
% %plot based on pulses
% figure;
% p = panel;
% p.pack(1,2)
% p(1,1).select();
% hold on
% for AMP=1:length(AMPall)
%     errorbar(Pulseall,ratepeak_suppression(AMP,:),ratepeakstd_suppression(AMP,:),'Color',[0 AMP/length(AMPall) 1/AMP])
% end
% lgd=legend('2','5','6','8','10');
% lgd.Title.String = 'Current \muA';
% title('Select data NOT stitch suppression Peak response')
% ylabel('Sp/s')
% xlabel('# Pulses')
% set(gca,'TickDir','out');
% %for p2 plot the average
% p(1,2).select();
% hold on
% for AMP=1:length(AMPall)
%     errorbar(Pulseall,ratespiking_suppression_mean(AMP,:),ratespiking_suppression_std(AMP,:),'Color',[0 AMP/length(AMPall) 1/AMP])
% end
% lgd=legend('2','5','6','8','10');
% lgd.Title.String = 'Current \muA';
% title('Select data NOT stitch suppression Average 1-30ms after stim')
% ylabel('Sp/s')
% xlabel('# Pulses')
% set(gca,'TickDir','out');
% 
% end