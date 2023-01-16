function [stimChn,CHN]=loadstimulationchannels
tparams = dir('*_exp_datafile_*.mat');
if isempty(tparams)
    TrialParams = [];
    return
end
tparams = tparams.name;
TrialParams=load(tparams,'TrialParams');
TrialParams=TrialParams.TrialParams;
stimChn = unique(cell2mat(TrialParams(2:end,3)));
stimChn=stimChn(stimChn~=0);
load(tparams,'CHN');
CHN=CHN;
end