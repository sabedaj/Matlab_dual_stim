%% Initial Setup
FS = 30000;
filepath = pwd;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir([filepath, filesep, 'info.rhs']);
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
name = pwd;
[FP,name,ext] = fileparts(name);
name = name(1:end-14);
%% Deal with the digital lines
disp('Cleaning the digital lines');
fileinfo = dir('digitalin.dat');
nSam = fileinfo.bytes/2;
digin_fid = fopen('digitalin.dat','r');
digital_in = fread(digin_fid, nSam, 'uint16');
fclose(digin_fid);
stimDig = flip(find(digital_in == 1)); % Fix for finding the falling line instead of the rising line
visDig = flip(find(digital_in == 2)); % Fix for finding the falling line instead of the rising line
dt = diff(stimDig);
kill = dt == -1;
stimDig(kill) = [];
dt = diff(visDig);
kill = dt == -1;
visDig(kill) = [];
nStamps = max([length(stimDig); length(visDig)]);
time_stamps = nan(1,nStamps);
time_stamps(1,1:length(stimDig)) = flip(stimDig);
%time_stamps(2,1:length(visDig)) = flip(visDig);
% if ~isempty(time_stamps)
%     if isnan(time_stamps(2,2))
%         time_stamps = time_stamps(1,:);
%     else
%         time_stamps = time_stamps(2,:);
%     end
% end
% Correct for 200 us delay
loadDelay;
time_stamps = time_stamps + delay*(1e-6)*FS;
% At this point, we need to correct for jitter
d = Depth(E_Mapnumber); SKIP = 1;
try 
    loadStimChn;
catch 
    stimChn = 16;
end
trig = time_stamps;
if (SKIP == 0)
warning('off','signal:findpeaks:largeMinPeakHeight');
for t = 1:length(trig)
    v = loadRAW(trig(t));
    [~,adj] = findpeaks(abs(diff(v(d(stimChn(1)),:))),'MinPeakHeight',500);
    if isempty(adj)
        continue;
    end
    adj = adj(1) - (500*FS/1e3);
    trig(t) = trig(t) + adj;
    v = loadRAW(trig(t));
end
end
% Adjust trigger lines to match datafile
% loadNTrials; loadNREP;
% if length(trig) > n_Trials
%     % Delete trigger lines
%     TrialParams = loadTrialParams; AMP = loadAMP; DUR = loadDUR;
%     tParams = cell2mat(TrialParams(1:end,2));
%     for nS = 1:length(CHN)
%         for nD = 1:length(DUR)
%             for nA = 1:length(AMP)
%                 id = nA + (nD-1)*length(AMP) + (nS-1)*length(DUR)*length(AMP);
%                 tmp = find(tParams == id);
%                 mk = tmp(n_REP);
%                 trig(tmp(tmp > mk)) = NaN;
%             end
%         end
%     end
%     trig(isnan(trig)) = [];
% end
if trig(1)<30000
    trig(1)=[];
end
fprintf('There are %.0f digital lines here.\n',length(trig));
TP=loadTrialParams;
if length(trig)<(length(TP)/2)
    %Dealing with missing triggers- find first ‘no stim’ trial with a pulse
    %and then find next ‘no stim’ trial along and assign -500 to all trigs
    %in between, then shift by one value to add the missing trigger. The
    %-500 is used as a flag later on in the code.
    file = pwd;
    mDIR = dir('amplifier.dat');
    mNAME = mDIR.name;
    mFID = fopen([file filesep mNAME], 'r');
    fprintf('Correcting for missing triggers.\n');
    loadStimParams;
    StimParams(1,:) = [];
    StimParams(1:2:end,:) = [];
    NoStimTrials=find(cell2mat(StimParams(:,16))==-1);
    for i=length(NoStimTrials):-1:1
        if NoStimTrials(i)>length(trig)
            NoStimTrials(i)=[];
        end
    end
    Trigs_to_Check = trig(NoStimTrials);
    BIN = [-50 50];
    for t=1:length(Trigs_to_Check)
        OFFSET = cast(nChn*2*(Trigs_to_Check(t)+(BIN(1)*FS/1e3)),'int64');
        fseek(mFID,OFFSET,'bof'); % used to point to specific point in the data indicated by offset from begginning of file
        v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
        if any(v(5,:)>2000) || any(v(15,:)>2000) 
            Start_trig=find(trig==Trigs_to_Check(t-1));
            End_trig=find(trig==Trigs_to_Check(t));
            trig(Start_trig:End_trig)=-500;
            Insert_Trig=-500;
            trig = [trig(1,1:Start_trig) Insert_Trig trig(1,Start_trig+1:end)];
            Trigs_to_Check = trig(NoStimTrials);
        end
    end

end
if length(trig)~=size(TP,1)/2
    trig=[trig -500*ones(size(TP,1)/2-length(trig),1)'];
end
trig_fid = fopen([name '.trig.dat'],'w');
fwrite(trig_fid,trig,'double');
fclose(trig_fid);