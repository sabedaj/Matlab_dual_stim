function [spk_all,rate_all,Peak_latency_all]=SigmoidGenerator(chnnum,Spk_array,timestart,timeend)
trig=loadTrig(0);
avgtimebs=2;%10x longer in baseline 
[Spike_trialstruct,baslinespike_trialstruct] = OnlinesortTrials(trig,Spk_array,chnnum,timestart,timeend,avgtimebs);
loadAMP_all;
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
chn=unique(trialinfo(:,2));
chn(chn==0)=[];
count=0;
allsingletrialIDs=[];
IDs=cell(length(chn),1);
                    

for trialloop=1:numsim_elect:length(trialinfo(:,2))
    if sum(trialinfo(trialloop:trialloop+numsim_elect-1,2)==0)==numsim_elect-1
        count=count+1;
        for chnloop=1:length(chn)
            if any(trialinfo(trialloop:trialloop+numsim_elect-1,2)==chn(chnloop))
                IDs{chnloop}(count,1:2)=[trialinfo(trialloop,17) trialloop];
            end
        end
        %allsingletrialIDs(count)=trialloop;%index   trialinfo(trialloop,1);
    end
end
        tparams = dir('*_exp_datafile_*.mat');
        tparams = tparams.name;
        load(tparams,'simultaneous_stim');
spk_all=[];
for recordchn=1:length(chnnum)
    for stimchn=1:length(chn)
        IDs{stimchn}=sortrows(IDs{stimchn});
        IDs{stimchn}=IDs{stimchn}(any(IDs{stimchn},2),:);

        for ID=1:size(IDs{stimchn},1)
            t=['ID_' num2str((IDs{stimchn}(ID,2)))];
            spkcount=zeros(40,1);
            basespikecount=zeros(40,1);
            for j=1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn))]).(t),'AsArray',true),2)
                t2=['Trial_' num2str(j)];
                spkcount(j)=size(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn))]).(t).(t2),1);
                basespikecount(j)=size(baslinespike_trialstruct.(['Chn_' num2str(chnnum(recordchn))]).(t).(t2),1)/avgtimebs;
            end
           % p=ranksum(basespikecount,spkcount,'tail','left');
            basesubtractspike=spkcount-basespikecount;
            %if p<0.05
                spk_all.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))])(ID)=mean(basesubtractspike(1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn))]).(t),'AsArray',true),2)))/(timeend-timestart)*1000;
                [rate,Peak_latency]=spikelatency(chnnum(recordchn),Spk_array,IDs{stimchn}(ID,1:2),trig,chn(stimchn));
                rate_all.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))])(ID,:)=rate.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))]);
                Peak_latency_all.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))])(1,ID)=Peak_latency.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))]);
%             else
%                 spk_all.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))])(ID)=nan(1,1);
%                 rate_all.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))])(ID,:)=nan(1,181);
%                 Peak_latency_all.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))])(1,ID)=nan(1,1);
%             end
%             if ID==size(IDs{stimchn},1)
%                 p=ranksum(basespikecount,spkcount,'tail','left');
%                 if p>0.05
%                     spk_all.(['Chn_' num2str(chnnum(recordchn))]).(['stimchn_' num2str(chn(stimchn))])(1:size(IDs{stimchn},1))=nan(1,size(IDs{stimchn},1));
%                 end
%             end
        end
    end
end


end