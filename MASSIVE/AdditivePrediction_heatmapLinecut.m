function [predictdatastruct,ratio_all]=AdditivePrediction_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT, stderrspktrial,startpointseconds, secondstoanalyse, depthdriven)
%Creates a prediction of dual electrode stimulation based on single
%stimulation data to remove bias 

trialinfo=loadTrialInfo(0);
loadNORECORDELECT;
loadStimChn;
loadAMP_all;
AMPorig=loadAMP;
predictdatastruct=[];
AMP=AMP_all';
TrialParams=loadTrialParams;
[~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
[PlatChnArray, AmpPlatVal]=AllElectrodeResponseCurve_SM(AMPInterestSingleLinePlotINDEXDUAL,avgnospT,stderrspktrial,0);
AmpPlatVal=(1000/(secondstoanalyse-startpointseconds)).*AmpPlatVal;
avgnostim=mean(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1),2);%average reponse without stimulation
%find closest AMP level to your chosen amp level to analyse (in case you
%pick a level not tested)
AMPInterestSingleLinePlot=AMPInterestSingleLinePlot/2;
[~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will te ll you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
    originalEND=size(TrialParams,1);
else
    E_Mapnumber=loadMapNum;
    loadoriginalEND;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end

AMP_double=[0 AMP(2:end).*2];
AMP_quater=[0 AMP(2:end).*(2/4)];
AMP_threequaters=[0 AMP(2:end).*(2*3/4)];
Desired_trialthreequatersone=zeros(1,length(AMP));
Desired_trialquaterone=zeros(1,length(AMP));
Desired_trialthreequaterstwo=zeros(1,length(AMP));
Desired_trialquatertwo=zeros(1,length(AMP));
Desired_trialequal=zeros(1,length(AMP));
Desired_trialequaltwo=zeros(1,length(AMP));
for Chosenstimchn=1:length(CHN) %used for stimulating channel one
    desiredchanneltrial=find(cell2mat(trialinfo(:,2))==stimChn(Chosenstimchn)); %finds trials with desired initial electrode
    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds trials with matching recording electrode spacing
    Desired_trialfirst=cell2mat(trialinfo(desiredchannel__singleampmath,1));%array of mathcning trial number
    normalisedspkthreequaterschnone=zeros(nChn,length(Desired_trialfirst));
    normalisedspkquaterchnone=zeros(nChn,length(Desired_trialfirst));
    normalisedAvgspikingTchnone=zeros(nChn,length(Desired_trialfirst));
    for Amploop=1:length(AMP)%used to identify existing 75/25 amplitudes for chn 1
        try
            Desired_trialthreequatersone(Amploop)=(desiredchannel__singleampmath(cell2mat(trialinfo(desiredchannel__singleampmath,18))==AMP_threequaters(Amploop))+1)/2;%array of mathcning trial number
        catch
            Desired_trialthreequatersone(Amploop)=0;
        end
        try
            Desired_trialquaterone(Amploop)=(desiredchannel__singleampmath(cell2mat(trialinfo(desiredchannel__singleampmath,18))==AMP_quater(Amploop))+1)/2;%array of mathcning trial number
        catch
            Desired_trialquaterone(Amploop)=0;
        end
        try
            Desired_trialequal(Amploop)=(desiredchannel__singleampmath(cell2mat(trialinfo(desiredchannel__singleampmath,18))==AMP(Amploop))+1)/2;%array of mathcning trial number
        catch
            Desired_trialequal(Amploop)=0;
        end

    end
    normalisedAvgspikingTchnone(:,Desired_trialequal~=0)=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,Desired_trialequal(Desired_trialequal~=0))-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
    normalisedspkthreequaterschnone(:,Desired_trialthreequatersone~=0)=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,Desired_trialthreequatersone(Desired_trialthreequatersone~=0))-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
    normalisedspkquaterchnone(:,Desired_trialquaterone~=0)=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,Desired_trialquaterone(Desired_trialquaterone~=0))-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
    

    for singlecond=1:length(NORECORDELECT) %used for stimulating channel two
        
        desiredchanneltwotrial=find(cell2mat(trialinfo(:,2))==(stimChn(Chosenstimchn)+NORECORDELECT(singlecond)+1)); %finds trials with desired initial electrode
        desiredchanneltwo__singleampmath=desiredchanneltwotrial((cell2mat(trialinfo(desiredchanneltwotrial+1,2))==(0)),1);% finds trials with matching recording electrode spacing
        Desired_trialsecond=cell2mat(trialinfo(desiredchanneltwo__singleampmath,1));%array of mathcning trial number
        normalisedAvgspikingTchntwo=zeros(nChn,length(Desired_trialsecond));
        normalisedspkthreequaterschntwo=zeros(nChn,length(Desired_trialsecond));
        normalisedspkquaterchntwo=zeros(nChn,length(Desired_trialsecond));
        for Amploop=1:length(AMP) %used to identify existing 75/25 amplitudes for chn 2
            try
                Desired_trialthreequaterstwo(Amploop)=(desiredchanneltwo__singleampmath(cell2mat(trialinfo(desiredchanneltwo__singleampmath,18))==AMP_threequaters(Amploop))+1)/2;%array of mathcning trial number
            catch
                Desired_trialthreequaterstwo(Amploop)=0;
            end
            try
                Desired_trialquatertwo(Amploop)=(desiredchanneltwo__singleampmath(cell2mat(trialinfo(desiredchanneltwo__singleampmath,18))==AMP_quater(Amploop))+1)/2;%array of mathcning trial number
            catch
                Desired_trialquatertwo(Amploop)=0;
            end
            try
                Desired_trialequaltwo(Amploop)=(desiredchanneltwo__singleampmath(cell2mat(trialinfo(desiredchanneltwo__singleampmath,18))==AMP(Amploop))+1)/2;%array of mathcning trial number
            catch
                Desired_trialequaltwo(Amploop)=0;
            end
        end
        normalisedAvgspikingTchntwo(:,Desired_trialequaltwo~=0)=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,Desired_trialequaltwo(Desired_trialequaltwo~=0))-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
        normalisedspkthreequaterschntwo(:,Desired_trialthreequaterstwo~=0)=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,Desired_trialthreequaterstwo(Desired_trialthreequaterstwo~=0))-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
        normalisedspkquaterchntwo(:,Desired_trialquatertwo~=0)=(1000/(secondstoanalyse-startpointseconds))*(avgnospT(:,Desired_trialquatertwo(Desired_trialquatertwo~=0))-avgnostim); %normalises data by subtracting no stim trials and converts to spikes per second
      
        
        [bfequal]=regression3Dplot(avgnospT,stderrspktrial,startpointseconds, secondstoanalyse, AMPInterestSingleLinePlot, AMPInterestSingleLinePlot);
        avgratioE1E2equal=(((bfequal)/(bfequal+1))*2);%multiply by two so the (1-b) and B ratio works. divide by two for the average
        avgratioE1E2=1;
        Additivespkequal=(avgratioE1E2).*normalisedAvgspikingTchnone+(2-avgratioE1E2).*normalisedAvgspikingTchntwo;
        [bfthreequaters]=regression3Dplot(avgnospT,stderrspktrial,startpointseconds, secondstoanalyse, AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL),AMP_quater(AMPInterestSingleLinePlotINDEXDUAL));
        
        avgratioE1E2threequaters=(((bfthreequaters)/(bfthreequaters+1))*2);%multiply by two so the (1-b) and B ratio works. divide by two for the average
        avgratioE1E2=1;
        Additivespkthreequaters_Eone=(avgratioE1E2).*normalisedspkthreequaterschnone+(2-avgratioE1E2).*normalisedspkquaterchntwo;
        Additivespkthreequaters_Eone=Additivespkthreequaters_Eone(:,(Desired_trialquatertwo~=0)+(Desired_trialthreequatersone~=0)==2);%only includes columns that had additive values for both conditions
        Additivespkthreequaters_Eone=[zeros(nChn,1) Additivespkthreequaters_Eone];
        [bfquater]=regression3Dplot(avgnospT,stderrspktrial,startpointseconds, secondstoanalyse, AMP_quater(AMPInterestSingleLinePlotINDEXDUAL),AMP_threequaters(AMPInterestSingleLinePlotINDEXDUAL));
        avgratioE1E2quater=(((bfquater)/(bfquater+1))*2);%multiply by two so the (1-b) and B ratio works. divide by two for the average
        avgratioE1E2=1;
        Additivespkquater_Eone=(avgratioE1E2).*normalisedspkquaterchnone+(2-avgratioE1E2).*normalisedspkthreequaterschntwo;
        Additivespkquater_Eone=Additivespkquater_Eone(:,(Desired_trialthreequaterstwo~=0)+(Desired_trialquaterone~=0)==2);%only includes columns that had additive values for both conditions
        Additivespkquater_Eone=[zeros(nChn,1) Additivespkquater_Eone];
        
        ratio_all=[bfquater,bfequal,bfthreequaters];

        %used to check if the additive spiking level is above the plateu
        %for each channel
        for chncount=1:nChn
            for replaceval=1:size(Additivespkequal,2)
                if Additivespkequal(chncount,replaceval)>AmpPlatVal(chncount)+0.15*AmpPlatVal(chncount)
                    Additivespkequal(chncount,replaceval)=AmpPlatVal(chncount)+0.15*AmpPlatVal(chncount);
                end
            end
        end
        for chncount=1:nChn
            for replaceval=1:size(Additivespkthreequaters_Eone,2)
                if Additivespkthreequaters_Eone(chncount,replaceval)>AmpPlatVal(chncount)+0.15*AmpPlatVal(chncount)
                    Additivespkthreequaters_Eone(chncount,replaceval)=AmpPlatVal(chncount)+0.15*AmpPlatVal(chncount);
                end
            end
        end
        for chncount=1:nChn
            for replaceval=1:size(Additivespkquater_Eone,2)
                if Additivespkquater_Eone(chncount,replaceval)>AmpPlatVal(chncount)+0.15*AmpPlatVal(chncount)
                    Additivespkquater_Eone(chncount,replaceval)=AmpPlatVal(chncount)+0.15*AmpPlatVal(chncount);
                end
            end
        end
        check=['T50_50_' num2str(CHN(Chosenstimchn))];
        predictdatastruct.(check) = Additivespkequal;
        check=['T75_25_' num2str(CHN(Chosenstimchn))];
        predictdatastruct.(check) = Additivespkthreequaters_Eone;
        check=['T25_75_' num2str(CHN(Chosenstimchn))];
        predictdatastruct.(check) = Additivespkquater_Eone;
        %%%%PLOTTING
        
        %HEATMAPS
        %75E2/25E1
        chosen_trials=(Desired_trialfirst(2)+2)*ones(size(Additivespkquater_Eone,2),1);
        figure
        fignum=gcf;
        fignum=fignum.Number;
        if size(Additivespkquater_Eone,2)>1
            subplot(1,3,1)
            chosen_trials_amp=[0 AMP((Desired_trialthreequatersone~=0)+(Desired_trialquatertwo~=0)==2)];
            DepthChangeingSpiking_SM(Additivespkquater_Eone, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
            title(['ESTIMATION Stimchn: ' num2str(cell2mat(trialinfo((Desired_trialfirst(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((Desired_trialsecond(2)*2-1),2))) ' @ 25/75' ])
            xt = xticks;
            %labelxx=AMP((Desired_trialthreequatersone~=0)+(Desired_trialquatertwo~=0)==2);
            xticklabels(cellstr(string((xt.*2))));
            %xticklabels(cellstr(string((AMP_double((Desired_trialthreequatersone~=0)+(Desired_trialquatertwo~=0)==2)))));
        end
        
        %50/50
        chosen_trials=(Desired_trialfirst(2)+3)*ones(size(Additivespkequal,2),1);
        chosen_trials_amp=AMP((Desired_trialequal~=0)+(Desired_trialequaltwo~=0)==2);
        subplot(1,3,2)
        DepthChangeingSpiking_SM(Additivespkequal, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ESTIMATION Stimchn: ' num2str(cell2mat(trialinfo((Desired_trialfirst(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((Desired_trialsecond(2)*2-1),2))) ' @ 50/50' ])
        xt = xticks;
        xticklabels(cellstr(string((xt.*2))));
        xlim([0 AMP(end)/2]);
        %25E2/75E1
        chosen_trials=(Desired_trialfirst(2)+4)*ones(size(Additivespkthreequaters_Eone,2),1);
        if size(Additivespkthreequaters_Eone,2)>1 %if there is only one matching amplitude then you cannot generate a heatmap
            subplot(1,3,3)
            chosen_trials_amp=[0 AMP((Desired_trialthreequaterstwo~=0)+(Desired_trialquaterone~=0)==2)];
            DepthChangeingSpiking_SM(Additivespkthreequaters_Eone, chosen_trials, fignum,depthdriven, max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
            title(['ESTIMATION Stimchn: ' num2str(cell2mat(trialinfo((Desired_trialfirst(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((Desired_trialsecond(2)*2-1),2))) ' @ 75/25' ])
            xt = xticks;
            xticklabels(cellstr(string((xt.*2))));
            %xticklabels(cellstr(string((AMP_double((Desired_trialthreequaterstwo~=0)+(Desired_trialquaterone~=0)==2)))));
        end

        ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
        ax.Label.String='Sp/s';
        ax.Label.Rotation=270;
             
        
        %LINECUTS

        if (Desired_trialquatertwo(AMPInterestSingleLinePlotINDEXDUAL)~=0)&&(Desired_trialthreequatersone(AMPInterestSingleLinePlotINDEXDUAL)~=0)
            [~, AMPInterestSingleLinePlotINDEXDUALt]=min(abs(AMPorig(2:end)-AMPInterestSingleLinePlot*2)); %x is a throw away variable but it will tell you how close the value is to your chosen val
            AMPInterestSingleLinePlotINDEXDUALt=AMPInterestSingleLinePlotINDEXDUALt+1; %since you did not check -1 condition add one to the index AMP(2:end)
            Additivespkthreequaters_Eone=Additivespkthreequaters_Eone(:, AMPInterestSingleLinePlotINDEXDUALt);%only includes columns that had additive values for both conditions
        else 
            Additivespkthreequaters_Eone=zeros(nChn,1);
        end
        if (Desired_trialquaterone(AMPInterestSingleLinePlotINDEXDUAL)~=0)&&(Desired_trialthreequaterstwo(AMPInterestSingleLinePlotINDEXDUAL)~=0)
            [~, AMPInterestSingleLinePlotINDEXDUALt]=min(abs(AMPorig(2:end)-AMPInterestSingleLinePlot*2)); %x is a throw away variable but it will tell you how close the value is to your chosen val
            AMPInterestSingleLinePlotINDEXDUALt=AMPInterestSingleLinePlotINDEXDUALt+1; %since you did not check -1 condition add one to the index AMP(2:end)
             Additivespkquater_Eone= Additivespkquater_Eone(:, AMPInterestSingleLinePlotINDEXDUALt);%only includes columns that had additive values for both conditions
 
        else
            Additivespkquater_Eone=zeros(nChn,1);
        end
        Additivespkequal=Additivespkequal(:,AMPInterestSingleLinePlotINDEXDUAL);

        
        figure %plotting estimation line plot
        ax = gca;
        ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
        hold on
        SMOOTHING=1;
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        
        rate = conv(Additivespkquater_Eone,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
%rate = Additivespkquater_Eone;        
        plot(rate)

        rate = conv(Additivespkequal,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        %rate = Additivespkequal; 
        plot(rate)

        rate = conv(Additivespkthreequaters_Eone,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        %rate = Additivespkthreequaters_Eone;
        plot(rate)
        
        ylabel('Sp/s')
        xlabel('Channel number')
        xlim([1 nChn])
        xline((stimChn(Chosenstimchn)+NORECORDELECT(singlecond)+1),'-.k')
        xline((stimChn(Chosenstimchn)),'k')
        legend('75/25', '50/50','25/75', 'Stim E1','Stim E2')
        title([num2str(AMP_double(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((stimChn(Chosenstimchn)+NORECORDELECT(singlecond)+1)) ' & ' num2str(stimChn(Chosenstimchn)) ' ESTIMATION'])
    end
end

end

