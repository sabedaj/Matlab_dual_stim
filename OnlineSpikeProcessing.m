
%% Make sigmoid from chosen channels
%ONLY if there is single electrode stimulation
%can also plot raster - need to uncomment
electrodes='A-030,C-009';%"A-009'"% for all electrodes = Depth(1)
amplifier_channels=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);

% % all electrodes
% electrodes='';
% temp = ProbeMAP;
% for i=2:nChn+1
%     electrodes=strcat(electrodes,temp{i,6},',');
% end
% electrodes(end)=[];

Spk_array=OnlineSigmoidGenerator(electrodes,nChn);


%% firing rate response all channels 
trialinfo=loadTrialInfo;%use this to investigate which ID you want
ID=100;
amplifier_channels=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);
Chns='A-030,C-009';
Spk_array=OnlineFRresponse(Chns,nChn,ID);

%% heatmap
trialinfo=loadTrialInfo;%use this to investigate which ID you want
ID=5;
amplifier_channels=read_Intan_RHS2000_file;
nChn=length(amplifier_channels);
array=OnlineHeatmap(nChn,ID);