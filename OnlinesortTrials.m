function [Spike_trialstruct,baslinespike_trialstruct,latency_trialstruct] = OnlinesortTrials(trig,sp,chn,startpointms,mstoanalyse,varargin)
%Sort into trial IDs and plot example spikes
%OUTPUT - the structure containing spiking information for each trial/repeat where
%each cell relates to one trial ID. The array in this cell is in the format
%of (electrodes,trials/repeats) where the total number of spikes per
%repeat is calculated for the given time

%INPUT - starting point after trigegr to begin counting spikes
%(startpointms) in ms, end point after trigger to stop counting spikes
%(mstoanalyse) in ms, trigger information (trig)
FS=30000;
nstage=round(length(sp)./32);
nChn=length(sp);
TrialParams = loadTrialParams;
TrialParams=cell2mat(TrialParams);
maxtid=max(TrialParams(:,2));
nospI=[];
Spike_trialstruct=[];
baslinespike_trialstruct=[];
for headstage=0:nstage-1
    C1 = intersect(sp{1+(headstage*32)}(:,1),sp{2+(headstage*32)}(:,1));%checks later electrodes
    C2 = intersect(sp{24+(headstage*32)}(:,1),sp{6+(headstage*32)}(:,1)); %early electrodes
    C3 = intersect(sp{21+(headstage*32)}(:,1),sp{5+(headstage*32)}(:,1)); %middle electrodes
    test1= intersect(C1,C2);
    test2= intersect(C2,C3);
    check=intersect(test1,test2);
    for chncount=1:nChn
        [~,r,~]=(intersect(sp{chncount}(:,1),check));
        sp{chncount}(r,:)=[];
    end
end

% startpointms=2;
% mstoanalyse=8;
if nargin>5
    avgtimebs=varargin{1};
else
    avgtimebs=10;
end

% Sort into trial IDs
for tID=1:maxtid
    if isempty(TrialParams(TrialParams(1:end,2) == tID))
        fprintf('No trial ID: %d\n',tID)
    else
        TrialParamstID = find(TrialParams(1:end,2) == tID); %identifies trial row matching trial ID
        num_elect=min(diff(find(diff(TrialParamstID)~=1)));
        if num_elect>1
            TrialParamstID=TrialParamstID(num_elect:num_elect:end);
            TrialParamstID=TrialParamstID./num_elect;
        end
        trigtID = trig(TrialParamstID)./(FS/1000);
        trigtID(trigtID==-500/(FS/1000))=[];
        nTrig = length(trigtID);
        IDnum=['ID_' int2str(tID)];
        for indT=1:nTrig
            tnum=['Trial_' int2str(indT)];
            for chsp=1:1:length(chn)
                v = sp{chn(chsp)};
                if ~isempty(v)
                    spikedetailstrig=v(((v(:,1)>(trigtID(indT)+startpointms))&(v(:,1)<(trigtID(indT)+mstoanalyse))),:); %v(((v(:,1)>0)&(v(:,1)<20000)),:);
                    
                    %filter out noise
%                     indexdelete=false(size(spikedetailstrig,1),1);
%                     for i=1:size(spikedetailstrig,1)
%                         [~,indexspikemin]=min(spikedetailstrig(i,2:end));
%                         [val,indexspikemax]=max(spikedetailstrig(i,2:end));
%                         if indexspikemin~=13 || (indexspikemax>25 && val>150) || (spikedetailstrig(i,2)>75)
%                             indexdelete(i)=true;
%                         end
%                     end
%                     spikedetailstrig(indexdelete,:)=[];
%                     if any(spikedetailstrig(:,2:end)<-100,'all')
%                       figure(1);hold on; plot(1/30:1/30:49/30,spikedetailstrig(:,2:end))
%                     end

                    timems=mstoanalyse-startpointms;
                    baslinespiketrig=v(((v(:,1)>(trigtID(indT)-5-avgtimebs*timems))&(v(:,1)<(trigtID(indT)-5))),:); %v(((v(:,1)>0)&(v(:,1)<20000)),:);
                    latencymeasure=v(((v(:,1)>(trigtID(indT)-90))&(v(:,1)<(trigtID(indT)+90))),:) - trigtID(indT); %v(((v(:,1)>0)&(v(:,1)<20000)),:);
                    
                else
                    spikedetailstrig=[];
                    baslinespiketrig=[];
                    latencymeasure=nan(50,1);
                end
                mchsp=chn(chsp);
                cnum=['Chn_' int2str(mchsp)];
                Spike_trialstruct.(cnum).(IDnum).(tnum)=spikedetailstrig;
                baslinespike_trialstruct.(cnum).(IDnum).(tnum)=baslinespiketrig;
                latency_trialstruct.(cnum).(IDnum).(tnum)=latencymeasure(:,1);
            end
        end
        
    end
end

end

