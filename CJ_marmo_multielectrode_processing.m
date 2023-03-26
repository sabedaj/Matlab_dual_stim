%multielectrode processing
%% loop through data 
D_data=dir;
namedat={D_data.name};
for k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    try
        cd([D_data(k).folder filesep currD])
        filepath=pwd;
        [filepathm,name,ext] = fileparts(filepath);
        %amp(k,1:5)=loadAMP;
        sp=load([name(1:end-14) '.sp.mat'],'sp');
        sp=sp.sp;
    catch
        continue
    end
        [stimChn,~]=loadstimulationchannels;
    chn_range=1:128;
        if str2double(name(end-12:end-7))<220812
        warning('port D bad')
        chn_range(chn_range>96)=[];
        end
     timestart=2;
    timeend=90;% time to stop looking for spikes in ms
    trig=loadTrig(0);
              TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
       [IDstruct, baslinespikestruct] = sortTrials_SM(timestart,timeend,trig,0,1,1,endtrial);
    [avgnospT,stderrspktrial,~] = AverageTrialResponse_SM(IDstruct, baslinespikestruct);
       
       sortCJelectrodeOrder_V1V2;
end

% sort data according to channel - measure rate over 100ms
