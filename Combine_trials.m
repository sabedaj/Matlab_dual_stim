function Combine_trials(Subdir1, Subdir2, varargin)
%%need to combine sp so that we can see V1 vs V2 activation. also run the
%%visual stimuli

SubDir={Subdir1;Subdir2};
for loopdir=1:nargin-2
    SubDir{loopdir+2}=varargin{loopdir};
end
skip=1;
if skip~=1
concat_Trials=[];

for Trialdir=1:size(SubDir,1)
    cd(SubDir{Trialdir})
    load('IDstruct.mat','IDstruct','baslinespikestruct')
    if Trialdir==1
        concat_Trials=IDstruct;
        concat_Basline=baslinespikestruct;
    else
        for trial=1:length(fieldnames(IDstruct))
            tid=['T' num2str(trial)];
            concat_Trials.(tid)=[concat_Trials.(tid),IDstruct.(tid)];
            concat_Basline.(tid)=[concat_Basline.(tid),baslinespikestruct.(tid)];
        end
    end
end
save('Combined_trials.mat',"concat_Basline","concat_Trials")
[avgnospT,stderrspktrial,trialinfo]=AverageTrialResponse_SM(concat_Trials,concat_Basline);
end
%% raster combine

for Trialdir=1:size(SubDir,1)
    cd(SubDir{Trialdir})
    [chninfo,~]=read_Intan_RHS2000_file;
    sp=loadSpikes;
    nChn=length(chninfo);
    trig = loadTrig(0);
    trialinfo=loadTrialInfo;
    trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
    numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
    TP = loadTrialParams;
    for chn=1:nChn
        chnnum=['C' num2str(chn)];
        Spike_array=sp{chn};
        for ID=1:max(TP(:,2))
            tID = find(cell2mat(TP(numsim_elect:numsim_elect:end,2)) == ID);
            ID_trial=['T' num2str(ID)];
            theseTrig = trig(tID)./30;
            nT=length(theseTrig);
            %% Set up the raster data structure
            BIN = [-98 98];
           MAX = 400;
            xdata = [];
            ydata = [];
            for tr = 1:nT
                theseSp = (Spike_array(Spike_array(:,1) > theseTrig(tr)+BIN(1) & Spike_array(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
                for i = 1:length(theseSp)
                    xdata = [xdata, (theseSp(i))]; 
                    ydata = [ydata, tr*(MAX/nT)];
                end
            end
            XDATA_combine.(chnnum).(ID_trial)=xdata;
            YDATA_combine.(chnnum).(ID_trial)=ydata;
        end
    end
end
save('spikeinfo.mat',"YDATA_combine","XDATA_combine")
end