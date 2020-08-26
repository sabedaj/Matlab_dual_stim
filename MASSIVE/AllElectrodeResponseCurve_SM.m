function AllElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,h)
% Response curve at stimulating electrode

%INPUT - the trial information array, the average number of spikes per
%trial array, and the standard error per trial array. An array of ones and
%zeros denoting whether a channel had a significant increase in spiking

filepath = pwd;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir([filepath, filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=loadMapNum;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
E_MAP = Depth(E_Mapnumber);
loadStimChn;
AMP=loadAMP;
counter=0;

for count=1:length(CHN) %loop to find the correct trials
    desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)); %finds trials with desired initial electrode
    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
    Desired_trial=cell2mat(trialinfo(desiredchannel__singleampmath,1));
    counter=counter+1;
    if length(Desired_trial)>length(AMP)
        loopnum=length(Desired_trial)/length(AMP);
        temp_resultarray=zeros(loopnum,length(AMP));
        for reaarange=1:loopnum
            temp_resultarray(reaarange,1:length(AMP))=Desired_trial(1:length(AMP));
            Desired_trial(1:length(AMP))=[];
        end
        result_trial_array(counter:counter+1,1:length(AMP))=temp_resultarray;
        counter=counter+1;
    else
        result_trial_array(count,1:length(Desired_trial))=Desired_trial;%trials for each current Amp
    end
end

figure
if sum(h)==nChn
    hold on
    result_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    error_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    for Chncount=1:nChn%loop for plotting
        if h(Chncount)==1
            for i=1:size(result_trial_array,1)
                for j=1:size(result_trial_array,2)
                    result_array(i,j)=avgnospT((Chncount),result_trial_array(i,j));
                    error_array(i,j)=stderrspktrial((Chncount),result_trial_array(i,j));
                end
            end
            mean_result_array=mean(result_array,1);
            mean_error_array=mean(error_array,1);
            
            hold on
            errorbar(AMP,mean_result_array,mean_error_array)
            title('Significant Channel Comparison')
            ylabel('Average number of spikes')
            xlabel('Current uA')
        end
    end
    legend_names=1:nChn;
    legend(strsplit(num2str(legend_names)))
else
    subplot(1,2,1)
    hold on
    result_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    error_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    for Chncount=1:nChn%loop for plotting
        if h(Chncount)==1
            for i=1:size(result_trial_array,1)
                for j=1:size(result_trial_array,2)
                    result_array(i,j)=avgnospT((Chncount),result_trial_array(i,j));
                    error_array(i,j)=stderrspktrial((Chncount),result_trial_array(i,j));
                end
            end
            mean_result_array=mean(result_array,1);
            mean_error_array=mean(error_array,1);
            
            hold on
            errorbar(AMP,mean_result_array,mean_error_array)
            title('Significant Channel Comparison')
            ylabel('Average number of spikes')
            xlabel('Current uA')
        end
    end
    legend_names=1:nChn;
    legend_names=legend_names(logical(h));
    legend(strsplit(num2str(legend_names)))
    
    subplot(1,2,2)
    hold on
    result_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    error_array=zeros(size(result_trial_array,1),size(result_trial_array,2));
    for Chncount=1:nChn%loop for plotting
        if h(Chncount)==0
            for i=1:size(result_trial_array,1)
                for j=1:size(result_trial_array,2)
                    result_array(i,j)=avgnospT((Chncount),result_trial_array(i,j));
                    error_array(i,j)=stderrspktrial((Chncount),result_trial_array(i,j));
                end
            end
            mean_result_array=mean(result_array,1);
            mean_error_array=mean(error_array,1);
            
            hold on
            errorbar(AMP,mean_result_array,mean_error_array)
            title('Not Significant Channel Comparison')
            ylabel('Average number of spikes')
            xlabel('Current uA')
        end
    end
    legend_names=1:nChn;
    h=abs(h-1);
    legend_names=legend_names(logical(h));
    legend(strsplit(num2str(legend_names)))
end
end

