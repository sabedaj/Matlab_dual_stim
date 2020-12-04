function [normspk_all]=PoolNormalisedActivity(avgnospT,startpointseconds, secondstoanalyse)
%to calculate the peak of activity across all animals


trialinfo=loadTrialInfo(0);
lastwarn('', '');
loadNORECORDELECT;
loadStimChn;
loadVarAmp;
loadNREP;
TrialParams=loadTrialParams;
loadAMP_all;
AMP=AMP_all;
AMP_original=loadAMP;
normspk_all=[];
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

trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;%trials to jump over
endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions
if isempty(endtrialelect)
    endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
end
maxid=(originalEND-1)/(n_REP_true*2);
for AMPit=2:length(AMP_original)
    stimshankcentroid=zeros(20,1);
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
                %find closest AMP level to your chosen amp level to analyse (in case you
                %pick a level not tested)
                [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP_original(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
                AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
            else
                chosen_trials=Desired_trialequal;
                checkconsecutive=1;
                %find closest AMP level to your chosen amp level to analyse (in case you
                %pick a level not tested)
                [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
                AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
            end
            chosen_trials(chosen_trials>endtrialelect*(group_related))=[]; % removes any trials that are greater than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
            chosen_trials(chosen_trials<group_related/2)=[]; % removes any trials that are smaller than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
            normalisedAvgspikingT=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,chosen_trials));%-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
            loopcounter=loopcounter+1;
            singleLineplotval(:,loopcounter)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUAL);
        end
        SMOOTHING=1;
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        stackedacross=[];
        acrossshankplot=zeros(5,4);
        for shankplot=1:4
            
            for p=1:size(singleLineplotval,2)
                maxplot=max(singleLineplotval,[],2);
                maxplot(maxplot<0)=0;
                normspk=singleLineplotval(1+((shankplot-1)*16):(shankplot*16),p)./maxplot(1+((shankplot-1)*16):(shankplot*16));%varnostim(1+((shankplot-1)*16):(shankplot*16))';
                normspk(isnan(normspk))=0;
                normspk(isinf(normspk))=0;
                normspk(normspk<-1)=-1;
                rate = conv(normspk,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
                rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                acrossshankplot(p,shankplot)=mean(rate);
                stackedacross(1+((p-1)*16):(p*16),shankplot)=rate;%singleLineplotval(1+((shankplot-1)*16):(shankplot*16),p);
                rate1=rate;
                rate1(rate1<0)=0;
                A = trapz(1:16, rate1);
                B(1)=0;
                
                for lims = 2:16
                    B(lims) =  trapz(1:lims, rate1(1:lims));
                end
                
                [~,electrodecentroid]=min(abs(B-(A/2)));
                stimshankcentroid(p+(shankplot-1)*5)=electrodecentroid;
                check=['S' num2str(shankplot) '_T' num2str(p)];%not sorted
                if ~isempty(normspk_all)&&isfield(normspk_all,(check))
                    normspk_all.(check)=[normspk_all.(check),normspk];
                else
                    normspk_all.(check)=normspk;
                end
                
            end
        end
        
    end
    
end



end