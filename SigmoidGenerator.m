function [spk_all,rate_all,Peak_latency_all]=SigmoidGenerator(chnnum,Spk_array,timestart,timeend)
trig=loadTrig(0);
avgtimebs=10;%10x longer in baseline
[Spike_trialstruct,baslinespike_trialstruct,latency_trialstruct] = OnlinesortTrials(trig,Spk_array,chnnum,timestart,timeend,avgtimebs);
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
    chnname=['Chn_' int2str(chnnum(recordchn))];
    for stimchn=1:length(chn)
        stimchnname=['stimchn_' int2str(chn(stimchn))];
        IDs{stimchn}=sortrows(IDs{stimchn});
        IDs{stimchn}=IDs{stimchn}(any(IDs{stimchn},2),:);
        maxID=IDs{stimchn}(IDs{stimchn}(:,1)==10,2);
        p=0;
        checkbasesp=nan(40,length(IDs{stimchn}));
        checkspkcount=nan(40,length(IDs{stimchn}));
        for ID=1:size(IDs{stimchn},1)
            t=['ID_' int2str((IDs{stimchn}(ID,2)))];
            spkcount=zeros(40,1);
            basespikecount=zeros(40,1);
            timespikes=[];
            nT=size(struct2table(Spike_trialstruct.(chnname).(t),'AsArray',true),2);
            for j=1:nT
                t2=['Trial_' int2str(j)];
                spkcount(j)=size(Spike_trialstruct.(chnname).(t).(t2),1);
                basespikecount(j)=size(baslinespike_trialstruct.(chnname).(t).(t2),1)/avgtimebs;
                timespikes=[timespikes; latency_trialstruct.(chnname).(t).(t2)];
            end
            
            checkbasesp(1:40,ID)=basespikecount;
            checkspkcount(1:40,ID)=spkcount;
            if IDs{stimchn}(ID,2)==maxID
            p=ranksum(checkbasesp(:),checkspkcount(:),'tail','left');
            end
            basesubtractspike=spkcount-basespikecount;
%             if p>0.05
%                 spk_all.(chnname).(stimchnname)(1:size(IDs{stimchn},1))=nan(size(IDs{stimchn},1),1);
%                 rate_all.(chnname).(stimchnname)(1:size(IDs{stimchn},1),:)=nan(size(IDs{stimchn},1),181);
%                 Peak_latency_all.(chnname).(stimchnname)(1,1:size(IDs{stimchn},1))=nan(size(IDs{stimchn},1),1);
%                 continue
%             end
                
            spk_all.(chnname).(stimchnname)(ID)=mean(basesubtractspike(1:size(struct2table(Spike_trialstruct.(chnname).(t),'AsArray',true),2)))/(timeend-timestart)*1000;
            rate=hist(timespikes,-90:90);
            rateh=rate.*1000/nT;%rate.Values.*1000/nT;
            rate_all.(chnname).(stimchnname)(ID,:)=rateh;
            [~,peak]=max(rateh(90:end));
            Peak_latency_all.(chnname).(stimchnname)(1,ID)=peak;

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