function multipulseplotdata(ratestruct,savefilename,area)
s = RandStream('mt19937ar','Seed',296);
AMPall=[2 5 6 8 10];
baselinetime=1:89;
Pulseall=1:5;
% pre-processing
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=[];
if strcmp(area,'v1')
    chnrng=65:128;
elseif strcmp(area,'v2')
    chnrng=1:64;
end

for numerfolders=1:numfolderstotal
    trialend=size(savefilename{numerfolders}{4},1);
    for trial=1:trialend
        trialcheck=['T', num2str(trial)];
        
        AMPcheck=['A' num2str(savefilename{numerfolders}{4}(trial,2))];
        if strcmp(AMPcheck,'A-1')
            continue; %AMPcheck='A0';
        end
        Pulsecheck=['P' num2str(savefilename{numerfolders}{4}(trial,4))];
        if isfield(ratespiking,(AMPcheck)) && isfield(ratespiking.(AMPcheck),(Pulsecheck))
            ratespiking.(AMPcheck).(Pulsecheck)=cat(1,ratespiking.(AMPcheck).(Pulsecheck),nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,:,:).*1000,3));%-nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,baselinetime,:).*1000,[2 3]));
            flnm.(AMPcheck).(Pulsecheck)=[flnm.(AMPcheck).(Pulsecheck);savefilename{numerfolders}{5}];%-nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,baselinetime,:).*1000,[2 3]));
        else
            ratespiking.(AMPcheck).(Pulsecheck)=nanmean(ratestruct{numerfolders}.(trialcheck)(chnrng,:,:).*1000,3);
            flnm.(AMPcheck).(Pulsecheck)=[savefilename{numerfolders}{5}];
        end
        
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
        subselectbaseline=1:2:89;
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
        
        ratespiking_mean(AMP,Pulse)=nanmean(ratespiking.(AMPcheck).(Pulsecheck)(:,91:110),'all');
        ratespiking_std(AMP,Pulse)=SEM(ratespiking.(AMPcheck).(Pulsecheck)(:,91:110),0);
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
    title('Average 1-20ms after stim')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');
figure;
    hold on;
for Pulse=1:length(Pulseall)
    if Pulse>1
errorbar(AMPall,sum(ratespiking_mean(:,1:Pulse),2),ratespiking_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
    else
        errorbar(AMPall,(ratespiking_mean(:,Pulse)),ratespiking_std(:,Pulse),'Color',[0 Pulse/length(Pulseall) 1/Pulse])
    end
end
 lgd=legend('1','2','3','4','5');
lgd.Title.String = '# pulses';
    title('Average 1-20ms after stim,Pulses added')
    ylabel('Sp/s')
    xlabel('Current \muA')
    set(gca,'TickDir','out');

end