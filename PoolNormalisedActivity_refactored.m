function [Csplit_shankdist, Csplit_depthsep, Csplit]=PoolNormalisedActivity_refactored(avgnospT,startpointseconds, secondstoanalyse, ElectLayerClass)
%load info+data
trialinfo=loadTrialInfo(0);
loadStimChn;
if stimChn(1)<17 %determines the shank with the stimulating electrodes
    shank=1;
elseif stimChn(1)<33 && stimChn(1)>16
    shank=4;%4 but for the purpose of removing stim shank==2
    stimChn=stimChn-16;
elseif stimChn(1)<49 && stimChn(1)>32
    shank=2;%2
    stimChn=stimChn-32;
else
    shank=3;%3
    stimChn=stimChn-48;
end

loadVarAmp;
TrialParams=loadTrialParams;
loadAMP_all;
AMP=AMP_all;
AMP_original=loadAMP;

filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
    originalEND=size(TrialParams,1);
    AMP=AMP_original;
else
    E_Mapnumber=loadMapNum;
    loadoriginalEND;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
trialjump=find(cell2mat(trialinfo(:,2))==0)/2;%trials to jump over
trialjump=trialjump(2);%identify second zero corresponding with second channel single sitm
endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions
if isempty(endtrialelect)
    endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
end
maxid=(originalEND-1)/(n_REP_true*2);

AMP_zero=AMP_original;
AMP_zero(AMP_zero==-1)=0;
%setup data structures
for trial=1:5
    for current = 1:length(AMP_zero)
        Csplit.(['C' num2str(AMP_zero(current))]).(['T' num2str(trial)])=[];
        for shanksep=0:3
            Csplit_shankdist.(['C' num2str(AMP_zero(current))]).(['T' num2str(trial)]).(['D' num2str(shanksep)])=[];
            Csplit_depthsep.(['C' num2str(AMP_zero(current))]).(['T' num2str(trial)]).(['D' num2str(shanksep)])=[];
        end
    end
end


for AMPit=1:length(AMP_original)
    loopcounter=0;
    checkconsecutive=1;
    singleLineplotval=zeros(nChn,trialjump);
    AMPInterestSingleLinePlot=AMP_original(AMPit);

    for group_related=1:endtrialelect*2:maxid*2 %group_related is used to go through groups of related trials
        for TJ_related=1:2:trialjump*2 %goes through trials related by trial jump
            desiredchanneltrial_one=(find((cell2mat(trialinfo(:,2))==cell2mat(trialinfo(TJ_related+(group_related-1),2))))+1)/2; %finds trials with desired initial electrode
            desiredchanneltrial_two=find(cell2mat(trialinfo(:,2))==cell2mat(trialinfo(TJ_related+1+(group_related-1),2)))/2; %finds trials with desired second electrode
            chosen_trials=intersect(desiredchanneltrial_one,desiredchanneltrial_two); % finds the trials that intersect with two trials of interest
            Desired_trialequal=zeros(1,length(AMP));
            for Amploop=1:length(AMP) %used to identify existing 75/25 amplitudes for chn 2
                try
                    Desired_trialequal(Amploop)=(chosen_trials(cell2mat(trialinfo(chosen_trials*2,18))==AMP(Amploop)));%array of mathcning trial number
                catch
                    Desired_trialequal(Amploop)=0;
                end
            end
            if all(chosen_trials<=maxid)&&(VarAmp==1)&&(~isempty(chosen_trials(diff(chosen_trials(chosen_trials<maxid))==1))) %deals with varying amplitude having the same two electrodes for all middle trials
                chosen_trials=chosen_trials(checkconsecutive:3:length(chosen_trials));
                checkconsecutive=checkconsecutive+1; %flag for middle trials
                if AMPit==1
                    chosen_trials=chosen_trials(1);
                end
                %find closest AMP level to your chosen amp level to analyse (in case you
                %pick a level not tested)
                [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP_original(1:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
            else
                chosen_trials=Desired_trialequal;
                if AMPit==1 && AMP_original(1)==-1
                    chosen_trials=chosen_trials(1);
                end

                %find closest AMP level to your chosen amp level to analyse (in case you
                %pick a level not tested)
                [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(1:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
                checkconsecutive=1;
            end
            chosen_trials(chosen_trials>endtrialelect*(group_related))=[]; % removes any trials that are greater than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
            chosen_trials(chosen_trials<group_related/2)=[]; % removes any trials that are smaller than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
            normalisedAvgspikingT=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,chosen_trials));%-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
            loopcounter=loopcounter+1;
            singleLineplotval(:,loopcounter)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUAL);
        end
        for shankplot=1:4
            for trial=1:size(singleLineplotval,2)
                if shankplot==2
                    realshank=4;% sorted
                elseif shankplot==3
                    realshank=2;% sorted
                elseif shankplot==4
                    realshank=3;% sorted
                else
                    realshank=1;% sorted
                end
                AMP_noneg=AMPInterestSingleLinePlot;
                if AMPInterestSingleLinePlot==-1
                    AMP_noneg=0;
                end
                normspk=singleLineplotval(1+((shankplot-1)*16):(shankplot*16),trial);
                shanklayclass=ElectLayerClass(1+((shankplot-1)*16):(shankplot*16));
                Csplit.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)])=[Csplit.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]) [NaN(16-stimChn(1),1); normspk; NaN(32-(16-stimChn(1))-16,1)]];
                Shanksep=realshank-shank;
                Csplit_shankdist.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))])=[Csplit_shankdist.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))]) [NaN(16-stimChn(1),1); normspk; NaN(32-(16-stimChn(1))-16,1)]];
                WMborder=find(shanklayclass==4,1,'last');
                if isempty(WMborder)
                    GIborder=find(shanklayclass==3,1,'last');
                    if isempty(GIborder)
                        SGborder=find(shanklayclass==2,1,'last');
                        Csplit_depthsep.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))])=[Csplit_depthsep.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))]) [NaN(24-SGborder,1); normspk; NaN(32-(24-SGborder)-16,1)]];%pos 24 is gran pos 25 is supra
                    else
                        Csplit_depthsep.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))])=[Csplit_depthsep.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))]) [NaN(21-GIborder,1); normspk; NaN(32-(21-GIborder)-16,1)]];%pos 21 is infra, pos 22 is gran
                    end
                else 
                    Csplit_depthsep.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))])=[Csplit_depthsep.(['C' num2str(AMP_noneg)]).(['T' num2str(trial)]).(['D' num2str(abs(Shanksep))]) [NaN(7-WMborder,1); normspk; NaN(32-(7-WMborder)-16,1)]];%pos 7 is WM, pos 8 is infra
                end
            end
        end
    end
end
end