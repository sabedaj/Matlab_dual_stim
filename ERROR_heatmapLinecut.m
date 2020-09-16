function ERROR_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,depthdriven)
%Plots true data and prediction, then plots error
loadVarAmp;
AMP=loadAMP;
loadAMP_all;
loadStimChn;
TrialParams=loadTrialParams;
trialinfo=loadTrialInfo;
loadNORECORDELECT;
trialinfo(1,:)=[];
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

[~, AMPInterestSingleLinePlotINDEXDUAL]=min(abs(AMP(2:end)-AMPInterestSingleLinePlot)); %x is a throw away variable but it will tell you how close the value is to your chosen val
AMPInterestSingleLinePlotINDEXDUAL=AMPInterestSingleLinePlotINDEXDUAL+1; %since you did not check -1 condition add one to the index AMP(2:end)
truedatastruct=TrueData_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse,depthdriven);
predictdatastruct=AdditivePrediction_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT, stderrspktrial,startpointseconds, secondstoanalyse, depthdriven);
 trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;%trials to jump over
 endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions 
 if isempty(endtrialelect)
     endtrialelect=size(trialinfo,1)/2; %if there is only one set of conditions in the dataset
 end
 if endtrialelect>(originalTrialEND)
     endtrialelect=(originalTrialEND);
 end

ampindex50=(AMP==AMP_all.*2);
ampindex50=logical(sum(ampindex50,2));
ampindex50(1,1)=true;
figure
fignum=gcf;
fignum=fignum.Number;
for chosenstimchn=1:length(CHN)
    if VarAmp==1
        check=['T75_25_' num2str(CHN(chosenstimchn))];
        error75=truedatastruct.(check)-predictdatastruct.(check);
        check=['T25_75_'  num2str(CHN(chosenstimchn))];
        error25=truedatastruct.(check)-predictdatastruct.(check);
        check=['T50_50_'  num2str(CHN(chosenstimchn))];
        error50=truedatastruct.(check)-predictdatastruct.(check)(:,ampindex50');
        subplot(1,3,3)
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(4+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error75, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 75/25' ])
        caxis([-50 100])
        subplot(1,3,2)
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(3+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error50, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 50/50' ])
        caxis([-50 100])
        subplot(1,3,1)
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(2+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error25, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 25/75' ])
        caxis([-50 100])       
        ax=colorbar('Position',[0.93 0.1 0.03 0.85]);
        ax.Label.String='Sp/s';
        ax.Label.Rotation=270;
        figure %plotting estimation line plot
        ax = gca;
        ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
        hold on
        SMOOTHING=1;
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        
        rate = conv(error25(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        
        rate = conv(error50(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        
        rate = conv(error75(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        
        ylabel('Sp/s')
        xlabel('Channel number')
        xlim([1 nChn])
        xline((CHN(chosenstimchn)+NORECORDELECT(1)+1),'-.k')
        xline((CHN(chosenstimchn)),'k')
        legend('75/25', '50/50','25/75', 'Stim E1','Stim E2')
        title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' ERROR'])

    else
        check=['T50_50_' num2str(CHN(chosenstimchn))];
        error50=truedatastruct.(check)-predictdatastruct.(check)(:,ampindex50');
        figure
        fignum=gcf;
        fignum=fignum.Number;
        chosen_trials_amp=[0 AMP(2:end)];
        chosen_trials=(3+endtrialelect*(chosenstimchn-1)):trialjump:endtrialelect*chosenstimchn;
        DepthChangeingSpiking_SM(error50, chosen_trials, fignum,depthdriven,max(cell2mat(TrialParams(:,2))),chosen_trials_amp);%plots heatmap of estimated plots
        title(['ERROR Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo((chosen_trials(2)*2)-1,2))) ' ' num2str(cell2mat(trialinfo((chosen_trials(2)*2),2))) ' @ 50/50' ])
        caxis([-50 100])
        
         
        figure %plotting estimation line plot
        ax = gca;
        ax.ColorOrderIndex=ax.ColorOrderIndex+1;%skips single electrode colouring
        hold on
        SMOOTHING=1;
        window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
        rate = conv(error50(:,AMPInterestSingleLinePlotINDEXDUAL),window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
        rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        plot(rate)
        ylabel('Sp/s')
        xlabel('Channel number')
        xlim([1 nChn])
        xline((CHN(chosenstimchn)+NORECORDELECT(1)+1),'-.k')
        xline((CHN(chosenstimchn)),'k')

        legend('50/50','Stim E1','Stim E2')
        title([num2str(AMP(AMPInterestSingleLinePlotINDEXDUAL)) 'uA StimChn: ' num2str((CHN(chosenstimchn)+NORECORDELECT(1)+1)) ' & ' num2str(CHN(chosenstimchn)) ' ERROR'])

    end
end


end

