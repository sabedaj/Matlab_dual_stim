%% For quick analysis during experiment

% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-150; %removes spikes below this threshold
artefact_high=100; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=12; %How long after the trigger do you want to analyse spikes for(ms)? 
printspiking=0;
par=0;

%% 1. Blank stimulus

FS=30000;
filepath = pwd;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir([filepath,'\info.rhs']);
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
dName='amplifier';
vFID = fopen([filepath '\' dName '.dat'],'r');
mem = memory;
T = mem.MaxPossibleArrayBytes ./ (2 * nChn * 30000);
fileinfo = dir([filepath '\' dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * 30000);
if t_len < (Overall_time_to_analyse+Startpoint_analyse)
    error('Time to analyse exceeds the data recording period')
end
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256;
elseif Overall_time_to_analyse~=0
    T=Overall_time_to_analyse;
end
info = fileinfo.bytes/2;
nL = (ceil(info / (nChn*FS*double(T)))+1);
vblank=[];
BREAK = 1;
N=1;
denoiseIntan_sab(filepath, dName, T, par, Startpoint_analyse, Overall_time_to_analyse);
trig = loadTrig(0);
theseTrig = trig;

%% 2. Thresholds & Mu
allExtract_sab_1(dName,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');

%% Extract trial information
%trialinfo=TrialInfo_sab;

%% 4. Calculate Structure of sorted trials according to IDs
[IDstruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,Overall_time_to_analyse);

%% 5. Calculates template of trials and spiking responses (Output in true electrode order)
[avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct);
cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition 
avgnostim=(avgnospT(:,cell2mat(trialinfo(1:2:end,18))==-1));
testifsignificant=(avgnospT(:,cell2mat(trialinfo(1:2:end,18))~=-1));
checknumnovsstim=size(testifsignificant,2)/size(avgnostim,2);
[h,p]=ttest(avgnostim',testifsignificant(1:nChn,1:size(avgnostim,2))');
h(isnan(h))=0;

%% 6. Print all trial details
Legend=cell(nChn,1);
loadStimChn;
loadStimParams;
stimamplitudes=unique(cell2mat(StimParams(2:end,16)));
fprintf('The current amplitude in ascending order is: \n');
fprintf('%d\n',stimamplitudes);
fprintf('The stimulation channels were: \n');
fprintf('%d\n',stimChn);
%% Plots unfiltered data
tID=[36]; %trial ID of interest
Chan=[15];%channels of interest - recommend only one or two for raw
Chosen_trig=10;
Binstart=-250;
Binend=400;
stimConcatenate(tID,Chosen_trig,Chan,'RAW',Binstart,Binend);
title(['Chn ' ,num2str(Chan), '. Stimchn ' ,num2str(cell2mat(trialinfo(tID*2-1,2))), ' at ' ,num2str(cell2mat(trialinfo(tID*2,18))), 'uA raw'])
xlabel('Time (ms)')

%% Plots filtered data for selected channels and trial ID
Binstart=-250;
Binend=400;
tID=[30]; %trial ID of interest
Chan=[12];%channels of interest 
Chosen_trig=15;
stimConcatenate(tID,Chosen_trig,Chan,'MU',Binstart,Binend);
title(['Chn ' ,num2str(Chan), '. Stimchn ' ,num2str(cell2mat(trialinfo(tID*2-1,2))), ' at ' ,num2str(cell2mat(trialinfo(tID*2,18))), 'uA filtered'])
xlabel('Time (ms)')

%% 7. Plot stimulating electrode response curve
StimElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial);

%% 8. plotting average electrode response for all electrodes classes as significant and not significant
AllElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,h);

%% 9. plotting electrodes between stim and stim electrodes
InbetweenElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial);

%% Extracts responses from chosen single electrode stimulation for specific hannel
Chan=[15 16 17 18];%can input as many channels of interest as you want
stimChn=15;%only one stim chan
singleElectrodeResponseCurve_SM(trialinfo,avgnospT,stderrspktrial,Chan,stimChn);

%% Curve according to chosen trial with number of trials to skip
Chan=[15 16];%can input as many as you want
starttrial=3;%starting trial to analyse
trialjump=5;%number of trials to jump for next consecutive response
singleElectrodeResponseCurveAMPLITUDE_SM(trialinfo,avgnospT,stderrspktrial,Chan,starttrial,trialjump);

%% Plot response curve change over time (3D graph)
 TimeBegin=2; %How long after the trigger do you want skip spike analysis(ms)? 
 TimeEnd=20; %How long after the trigger do you want to stop spike analysis(ms)? 
 TimeStep=2; %Timestep of spike analysis between TimeBegin and TimeEnd(ms)? 
 starttrial=5;%starting trial to analyse
 electrode=16;
trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %number of trials to jump for next consecutive response equal to cond as default
 TimeChangeinSpiking_SM(TimeBegin, TimeEnd, TimeStep, electrode, starttrial, trialjump); %can input end trial but haven't tested
 %%
 starttrial=32;
 trialjump=find(diff(cell2mat(trialinfo(:,18))),1,'first')/2;
 endtrial=60;
DepthChangeingSpiking_SM(avgnospT, starttrial, trialjump, endtrial);

%% PCA spike sorting attempt
% CHNinterest=25;
% plotPCAspiking_sab(CHNinterest);


%% Deleting non-original files from filepath
%Delete_non_orig_sab;


%% Raster or histogram of firing rate
chn=13;
generate_StackedRaster_sab(chn);

%% MUA response

MUA_response_SM;
Binstart=-250;
Binend=400;
tID=[10]; %trial ID of interest
Chan=[20 21 22 23];%channels of interest
Chosen_trig=1;
stimConcatenate(tID,Chosen_trig,Chan,'MUA',Binstart,Binend);
title(['Chn ', num2str(Chan),' at ' string(trialinfo(tID*2,18)) 'uA MUA (filtered)'])
xlabel('Time (ms)')