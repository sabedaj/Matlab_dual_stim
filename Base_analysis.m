function Base_analysis(SubDir_Path)
%%for analysing data on massive and saving denoised files
tic
% Parameters to alter
Startpoint_analyse=0; %set to 0 for no input
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
artefact=-500; %removes spikes below this threshold
artefact_high=500; %removes spikes above this threshold
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)?
secondstoanalyse=12; %How long after the trigger do you want to analyse spikes for(ms)?
printspiking=0;
par=0;

%% 1. Blank stimulus

cd(SubDir_Path)%change working directory to where data is stored - currently manually input
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
filepath = pwd;
fileinfo = dir([filepath filesep 'info.rhs']);
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;
dName='amplifier';
mem_check=dir('amplifier.dat');
T = mem_check.bytes ./ (2 * nChn * FS);
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * FS);
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

denoiseIntan_sab(filepath, dName, T, par, Startpoint_analyse, Overall_time_to_analyse,FS);
cleanTrig_sabquick;
trig = loadTrig(0);
%% 2. Thresholds & Mu
allExtract_sab_1(dName,filepath,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');

%% 4. Calculate Structure of sorted trials according to IDs
starttrial=1;
trialjump=1;
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
endtrial=maxid;
[IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
save('IDstruct.mat', 'IDstruct','baslinespikestruct')

%% 5. Calculates template of trials and spiking responses (Output in true electrode order)
[avgnospT,stderrspktrial,trialinfo] = AverageTrialResponse_SM(IDstruct, baslinespikestruct);
save('Averagetrialresponse.mat','avgnospT','stderrspktrial')
toc
cd(folder)
fprintf(['End of Analysis for: ' SubDir_Path newline])
end

