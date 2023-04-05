function data = simpleDetrend(data,N,T,trig)
TrialParams = loadTrialParams;
numelect=find(diff(cell2mat(TrialParams(:,1)))~=0,1,'first');

TrialParams = cell2mat(TrialParams(1:numelect:end,2))';
nChn = size(data,1); FS = 30000;

    BIN_b = [-0.095*FS 0.095*FS]; % Amount to detrend in samples
    missed_trials = TrialParams(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
    missed_trig = trig(trig > ((N-1)*T*FS)-3000 & trig < ((N-1)*T*FS)+3000);
    missed_trials((missed_trig==-500))=[];
    missed_trig((missed_trig==-500))=[];
    
    trials = TrialParams(trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000);
    trig = trig(trig >= ((N-1)*T*FS)+3000 & trig <= (N*T*FS)-3000);
    if ~isempty(missed_trig)
        trials = [missed_trials trials];
        trig = [missed_trig trig];
    end
    if (N == 1)
        trig = trig - ((N-1)*T*FS);
    else
        trig = trig - ((N-1)*T*FS) + FS;
    end
    for t = 1:length(trig) 

        if trials(t) == 5
            pause(0.01);
        end
        thisTrig = trig(t);
        for c = 1:nChn

             % Detrend +-95ms - don't care about artefacts beyond that
            data(c,thisTrig+BIN_b(1):thisTrig+BIN_b(2)) = detrend(data(c,thisTrig+BIN_b(1):thisTrig+BIN_b(2)),2);

        end
    end

end