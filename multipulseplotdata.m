function multipulseplotdata(ratestruct,savefilename,area)
s = RandStream('mt19937ar','Seed',296);
AMPall=[2 5 6 8 10];
baselinetime=1:89;
Pulseall=1:5;
 subselectbaseline=1:2:89;
% pre-processing
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=[];
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
    for trial=1:trialend
        trialcheck=['T', num2str(trial)];
        
        AMPcheck=['A' num2str(savefilename{numerfolders}{4}(trial,2))];
        if strcmp(AMPcheck,'A-1')
            continue; %AMPcheck='A0';
        end
        Pulse=savefilename{numerfolders}{4}(trial,4);
        Pulsecheck=['P' num2str(Pulse)];
        
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
        datatoinclude=nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,:,:).*1000,3);
        if isfield(ratespiking,(AMPcheck)) && isfield(ratespiking.(AMPcheck),(Pulsecheck))
            ratespiking.(AMPcheck).(Pulsecheck)=cat(1,ratespiking.(AMPcheck).(Pulsecheck),datatoinclude);%-nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,baselinetime,:).*1000,[2 3]));
            flnm.(AMPcheck).(Pulsecheck)=[flnm.(AMPcheck).(Pulsecheck);savefilename{numerfolders}{5}];%-nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,baselinetime,:).*1000,[2 3]));
        else
            ratespiking.(AMPcheck).(Pulsecheck)=datatoinclude;
            flnm.(AMPcheck).(Pulsecheck)=[savefilename{numerfolders}{5}];
        end
         thresh=squeeze(mean(datatoinclude(:,subselectbaseline),2)+(std(datatoinclude(:,subselectbaseline),[],2).*2));
        notsigsingle=squeeze(mean(datatoinclude(:,92:181),2))<=thresh;
        Ratestitchspiking.(yrmnth).(AMPcheck).(chnstitch).(Pulsecheck)(notsigsingle,:)=nan;

    end
end
%to do
%1. plotting - check slack for details
%2. restrict to significant numbers only
ratespiking_mean=zeros(length(AMPall),length(Pulseall));
ratespiking_std=zeros(length(AMPall),length(Pulseall));
numsig=zeros(length(AMPall),5);
numelect=zeros(length(AMPall),5);
StimElectcount=zeros(length(AMPall),5);
%loop through amp and number of pulses then plot
avgtime=91:130;
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
       
        if isfield(ratespiking.(AMPcheck),(Pulsecheck))
        thresh=squeeze(mean(ratespiking.(AMPcheck).(Pulsecheck)(:,subselectbaseline),2)+(std(ratespiking.(AMPcheck).(Pulsecheck)(:,subselectbaseline),[],2).*2));
        notsigsingle=squeeze(mean(ratespiking.(AMPcheck).(Pulsecheck)(:,92:181),2))<=thresh;
        ratespiking.(AMPcheck).(Pulsecheck)(notsigsingle,:)=nan;
        numsig(AMP,Pulse)=sum(~notsigsingle);
        
        
        numelect(AMP,Pulse)= sum(~isnan(ratespiking.(AMPcheck).(Pulsecheck)(:,1)));
        % Reshape the array into a 2D matrix with 128 elements per column
        reshapedArray = reshape(ratespiking.(AMPcheck).(Pulsecheck)(:,1), 64, []);
        % Count the number of non-NaN values in each column
        StimElectcount(AMP,Pulse) = sum(~isnan(sum(~isnan(reshapedArray))));
        maxdat=max(ratespiking.(AMPcheck).(Pulsecheck)(:,avgtime),[],2);
        ratespiking_mean(AMP,Pulse)=nanmean(ratespiking.(AMPcheck).(Pulsecheck)(:,avgtime),'all');%nanmean(maxdat(maxdat~=0),'all');%nanmean(ratespiking.(AMPcheck).(Pulsecheck)(:,avgtime),'all');
        ratespiking_std(AMP,Pulse)=SEM(ratespiking.(AMPcheck).(Pulsecheck)(:,avgtime),0);
        stdshade(ratespiking.(AMPcheck).(Pulsecheck),0.2,[0 Pulse/length(Pulseall) 1/Pulse],[-90:300],1,ax)
        end
    end
    lgd=legend('1','2','3','4','5');
    lgd.Title.String = '# pulses';
    xlim([-90,90])
    title([num2str(AMPall(AMP)) '\muA'])
    ylabel('Sp/s')
    xlabel('Time (ms)')
    set(gca,'TickDir','out');
end

p(2,3).select();
hold on;
for Pulse=1:length(Pulseall)
errorbar(AMPall,ratespiking_mean(:,Pulse),ratespiking_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
end
 lgd=legend('1','2','3','4','5');
lgd.Title.String = '# pulses';
    title('Average 1-40ms after stim')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');

    %find nearest neighbour stim chn and add in the unique pulse data
    uniqueym=fields(Ratestitchspiking);
    for itYM=1:length(uniqueym)
        %iterate through Ratestitchspiking and average the data
        uniqueamp=fields(Ratestitchspiking.(uniqueym{itYM}));
        for itAMP=1:length(uniqueamp)
            uniquechn=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}));
            for itCHN=1:length(uniquechn)
                
                uniquepulse=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}));
                for itPulse=1:5
                    Pulsecheck=['P' num2str(itPulse)];
                    if ~isfield(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}),Pulsecheck) || size(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheck),2)~=391
                        %get number from uniquechn instead of 'CHNxxx'.
                        chnnum_all_text=regexp(uniquechn,'\d+','match');
                        % Convert all chnnums to numerical value
                        chnnum_all  = cellfun(@str2double, chnnum_all_text);
                        chninterest=chnnum_all(itCHN);
                        chnnum_all_text=uniquechn;
                        chnnum_all_text(itCHN)=[];
                        chnnum_all(itCHN)=[];
                        % find nearest neighbour stim chn to chninterest
                        for i=1:size(chnnum_all)
                            [~,idx] = min(abs(chninterest-chnnum_all));
                            if isfield(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}),Pulsecheck) && size(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck),2)==391
                                Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheck)=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck);
                                %add in pulse data from early pulses into late pulse
                                for pulsit=5:-1:itPulse+1
                                    Pulsecheckstitch=['P' num2str(pulsit)];
                                    if pulsit>1
                                    Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheckstitch)(:,94+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse))=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck)(:,94+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse));
                                    else
                                        Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(Pulsecheckstitch)(:,91+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse))=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(chnnum_all_text{idx}).(Pulsecheck)(chnrng,91+ceil(windowsize*(itPulse-1)):93+ceil(windowsize*itPulse));
                                    end
                                end
                                break;
                            else
                                chnnum_all_text{idx}=[];
                                chnnum_all(idx)=nan;
                                
                            end
                        end
                    end
                end
            end
        end
    end

    
    
    
    
    Ratestitchspikingmean=[];
    for itYM=1:length(uniqueym)
        %iterate through Ratestitchspiking and average the data
        uniqueamp=fields(Ratestitchspiking.(uniqueym{itYM}));
        for itAMP=1:length(uniqueamp)
            uniquechn=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}));
            for itCHN=1:length(uniquechn)
                uniquepulse=fields(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}));
                for itPulse=1:length(uniquepulse)
                    %put each unique combination of amp and pulse into a new array so that we can average them later ignoring uniqueym and chn
                    if length(uniquepulse)>3 && size(Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).P5,2)==391
                        if ~isfield(Ratestitchspikingmean,(uniqueamp{itAMP})) || ~isfield(Ratestitchspikingmean.(uniqueamp{itAMP}),(uniquepulse{itPulse}))
                            Ratestitchspikingmean.(uniqueamp{itAMP}).(uniquepulse{itPulse})=Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(uniquepulse{itPulse});
                        else
                            Ratestitchspikingmean.(uniqueamp{itAMP}).(uniquepulse{itPulse})=cat(1,Ratestitchspikingmean.(uniqueamp{itAMP}).(uniquepulse{itPulse}),Ratestitchspiking.(uniqueym{itYM}).(uniqueamp{itAMP}).(uniquechn{itCHN}).(uniquepulse{itPulse}));
                        end
                    end

                end
            end
        end
    
    end
    %plot Ratestitchspikingmean for each current with the lines on the plot being the pulse data
    figure;
    p = panel();
    p.pack(2, 3);
    ratestitchmean=zeros(length(AMPall),length(Pulseall));
    ratestitchstd=zeros(length(AMPall),length(Pulseall));
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
            if isfield(Ratestitchspikingmean,(AMPcheck)) && isfield(Ratestitchspikingmean.(AMPcheck),(Pulsecheck))
                stdshade(Ratestitchspikingmean.(AMPcheck).(Pulsecheck),0.2,[0 Pulse/length(Pulseall) 1/Pulse],[-90:300],1,ax)
                ratestitchmean(AMP,Pulse)=nanmean(Ratestitchspikingmean.(AMPcheck).(Pulsecheck)(:,avgtime),'all');
                ratestitchstd(AMP,Pulse)=SEM(Ratestitchspikingmean.(AMPcheck).(Pulsecheck)(:,avgtime),0);
            end
        end
        xlim([-90,90])
        title([num2str(AMPall(AMP)) '\muA'])
        ylabel('Sp/s')
        xlabel('Time (ms)')
        set(gca,'TickDir','out');
    end
    lgd=legend('1','2','3','4','5');
    lgd.Title.String = '# pulses';
p(2,3).select();
hold on
for Pulse=1:length(Pulseall)
    errorbar(AMPall,ratestitchmean(:,Pulse),ratestitchstd(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])

end
 lgd=legend('1','2','3','4','5');
lgd.Title.String = '# pulses';
    title('Average 1-20ms after stim,Pulses added')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');

end