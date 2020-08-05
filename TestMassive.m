% Parameters to alter
artefact=-100; %removes spikes below this threshold
artefact_high=100; %removes spikes above this threshold
startpointseconds=11; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=20; %How long after the trigger do you want to analyse spikes for(ms)? 
printspiking=0;
par=0;

%% 1. Blank stimulus
nChn=32;
FS=30000;
filepath = '/Volumes/shares/R-MNHS-Syncitium/Shared/Shared/Sabrina/Recordings/Rat_001/original/rat_001_005_005_200616_195625';
dName='amplifier';
vFID = fopen([filepath '\' dName '.dat'],'r');
mem = memory;
T = mem.MaxPossibleArrayBytes ./ (2 * 32 * 30000);
fileinfo = dir([filepath '\' dName '.dat']);
t_len = fileinfo.bytes/(32 * 2 * 30000);
if T > t_len
    T = t_len + 1;
end
if T > 256
    T = 256; 
end
info = fileinfo.bytes/2;
nL = (ceil(info / (nChn*FS*double(T)))+1);
vblank=[];
BREAK = 1;
N=1;
denoiseIntan_sab(filepath,dName,T,par)
trig = loadTrig(0);
theseTrig = trig;

%% 2. Thresholds & Mu
allExtract_sab_1(dName,T,par,artefact,artefact_high);% alternate -allExtract_sab(dName,T,par,artefact,artefact_high,trig,amp_issue);
fclose('all');