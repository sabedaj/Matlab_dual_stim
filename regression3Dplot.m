function [bh,bf]=regression3Dplot(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse)
%Creates Heatmaps of dual electrode stimulation and a linecut at the
%amplitude of interest input into the function

trialinfo=loadTrialInfo(0);
lastwarn('', '');
loadNORECORDELECT;
[~, warnId] = lastwarn();
if(~isempty(warnId))
    NORECORDELECT=[];
end


loadStimChn;
loadVarAmp;
loadNREP;
TrialParams=loadTrialParams;
loadAMP_all;
AMP=AMP_all;
AMP_original=loadAMP;
avgnostim=mean(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1),2);%average reponse without stimulation
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
truedatastruct=[];
 loopcounter=0;
 trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;%trials to jump over
 endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions 
 if isempty(endtrialelect)
     endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
 end
maxid=(originalEND-1)/(n_REP_true*2);
checkconsecutive=1;
singleLineplotval=zeros(nChn,trialjump);
singleLineplotvalHALVE=zeros(nChn,trialjump);
singleLineplotvalstd=zeros(nChn,trialjump);

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
            [~, halveAMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP_original(2:end)-AMPInterestSingleLinePlot/2)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             halveAMPInterestSingleLinePlotINDEXDUAL=halveAMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)

         else
             chosen_trials=Desired_trialequal;
             checkconsecutive=1;
             %find closest AMP level to your chosen amp level to analyse (in case you
             %pick a level not tested)
             [~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
            [~, halveAMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot/2)); %x is a throw away variable but it will tell you how close the value is to your chosen val
             halveAMPInterestSingleLinePlotINDEXDUAL=halveAMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
         end
         chosen_trials(chosen_trials>endtrialelect*(group_related))=[]; % removes any trials that are greater than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         chosen_trials(chosen_trials<group_related/2)=[]; % removes any trials that are smaller than the matching trial segment (e.g. stim chan 17 used as initial electrode in one set of trials and used as second electrode in second set trials)
         normalisedAvgspikingT=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,chosen_trials)-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
         stdsp=stderrspktrial(:,chosen_trials); %finds the standard deviation of chosen trials for plotting
         loopcounter=loopcounter+1;
         singleLineplotval(:,loopcounter)=normalisedAvgspikingT(:,AMPInterestSingleLinePlotINDEXDUAL);
         singleLineplotvalstd(:,loopcounter)=(1000/(secondstoanalyse-startpointseconds)).*stdsp(:,AMPInterestSingleLinePlotINDEXDUAL);

         singleLineplotvalHALVE(:,loopcounter)=normalisedAvgspikingT(:,halveAMPInterestSingleLinePlotINDEXDUAL);

         if cell2mat(trialinfo((chosen_trials(2)*2),2))==0 
             check=['T100_' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
             if TJ_related==1 && ~isempty(NORECORDELECT)

                 E1=num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)));
             elseif ~isempty(NORECORDELECT)
                E2=num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)));

             end
         else
             check=['T' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))*100/(cell2mat(trialinfo((chosen_trials(2)*2)-1,18))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2),18))*100/((cell2mat(trialinfo((chosen_trials(2)*2)-1,18)))+cell2mat(trialinfo((chosen_trials(2)*2),18)))) '_' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2)))];
         end
         %DepthChangeingSpiking_SM(normalisedAvgspikingT, chosen_trials,fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp); %plots heat maps
        
         
         truedatastruct.(check) = normalisedAvgspikingT;
         %          if VarAmp==1 && cell2mat(trialinfo(TJ_related+1+(TJ_related-1),2))==0
%              xlim([0 AMP(end)/2]) % makes axes comparable on single stim trials
%          end
     end
         figure %plotting real data line plot
         hold on
         scatter3(singleLineplotval(:,1),singleLineplotval(:,5),singleLineplotval(:,3))
         xlabel(['Channel ' E1 ' Sp/s'])
         ylabel(['Channel ' E2 ' Sp/s'])
         zlabel('E1+E2')
         title('Full')
         X_line=best_fit_3D_line([singleLineplotval(:,1),singleLineplotval(:,5),singleLineplotval(:,3)]);
         plot3(X_line(:,1),X_line(:,2),X_line(:,3),'-r','LineWidth',3) % best fit line 
         figure %plotting real data line plot
         hold on
         scatter3(singleLineplotvalHALVE(:,1),singleLineplotvalHALVE(:,5),singleLineplotval(:,3))
         X_line=best_fit_3D_line([singleLineplotvalHALVE(:,1),singleLineplotvalHALVE(:,5),singleLineplotval(:,3)]);
         plot3(X_line(:,1),X_line(:,2),X_line(:,3),'-r','LineWidth',3) % best fit line 
         xlabel(['Channel ' E1 ' Sp/s'])
         ylabel(['Channel ' E2 ' Sp/s'])
         zlabel('E1+E2')
         title('Half')
         
         figure 
         scatter(singleLineplotval(:,1),singleLineplotval(:,5))
         bf = singleLineplotval(:,1)\singleLineplotval(:,5);
         ycalc=bf*singleLineplotval(:,1);
         hold on
         plot(singleLineplotval(:,1),ycalc)
         xlabel(['Channel ' E1 ' Sp/s'])
         ylabel(['Channel ' E2 ' Sp/s'])
         title([num2str(AMPInterestSingleLinePlot) 'uA'])
         figure
         scatter(singleLineplotvalHALVE(:,1),singleLineplotvalHALVE(:,5))
         bh = singleLineplotvalHALVE(:,1)\singleLineplotvalHALVE(:,5);
         ycalc=bh*singleLineplotvalHALVE(:,1);
         hold on
         plot(singleLineplotvalHALVE(:,1),ycalc)
         xlabel(['Channel ' E1 ' Sp/s'])
         ylabel(['Channel ' E2 ' Sp/s'])
         title([num2str(AMPInterestSingleLinePlot/2) 'uA'])
         loopcounter=0;
end
end

