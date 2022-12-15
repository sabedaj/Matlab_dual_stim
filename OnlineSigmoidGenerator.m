function Spk_array=OnlineSigmoidGenerator(varargin)
%chns input 'A-001,B-013'
%nChn is the total number of channels recorded simulataneously
%% Online sigmoid generation
tic
filepath = pwd;
dName='amplifier';

%user input channels of choice?? 
if nargin==0
    chns=input('Please input channels to analyse like so "A-001,B-016":\n','s');
    nChn=input('Please input the total number of channels recorded from all headstages\n');
elseif nargin==1
    fprintf('Not enough input parameters\n')
    chns=input('Please input channels to analyse like so "A-001,B-016":\n','s');
    nChn=input('Please input the total number of channels recorded from all headstages\n');
else
    chns=varargin(1);
    nChn=varargin{2};
end
chnsplit=split(chns,",");
chnnum=zeros(size(chnsplit,1),1);
for i=1:size(chnsplit,1)
    if strcmp(chnsplit{i}(1),'A')
        chnnum(i)=str2double(chnsplit{i}(end-2:end))+1;
    elseif strcmp(chnsplit{i}(1),'B')
        chnnum(i)=str2double(chnsplit{i}(end-2:end))+33;
    elseif strcmp(chnsplit{i}(1),'C')
         chnnum(i)=str2double(chnsplit{i}(end-2:end))+65;
    elseif strcmp(chnsplit{i}(1),'D')
        chnnum(i)=str2double(chnsplit{i}(end-2:end))+97;
    end
end

%get small amplifier file
FS=30000;
vFID = fopen([filepath filesep dName '.dat'],'r');
BREAK=1;
fileinfo = dir([filepath filesep dName '.dat']);
t_len = fileinfo.bytes/(nChn * 2 * 30000);
T=256;
if t_len < T
   T=t_len;
end
temp=dir('*.trig.dat');
if isempty(temp)
    cleanTrig;
end
trig=loadTrig(0);
Spk_array=cell(size(chnnum,1),1,1);
N=0;
while BREAK
    v = fread(vFID,[nChn, (FS * T)],'int16') .* 0.195;
    N=N+1;
    if ~size(v,2)
        BREAK = 0;
    else
        trigN=trig(trig>T*FS*(N-1) & trig<=T*FS*N)-T*FS*(N-1);
        for chnit=1:length(chnnum)
            for trigit=1:length(trigN)
                if size(v,2)>=trigN(trigit)+FS*0.100
                    v(chnnum(chnit),trigN(trigit):trigN(trigit)+FS*0.100)=v(chnnum(chnit),trigN(trigit):trigN(trigit)+FS*0.100)-(mean(v(chnnum(chnit),trigN(trigit)+FS*0.005:trigN(trigit)+FS*0.025))-mean(v(chnnum(chnit),trigN(trigit)-FS*0.020:trigN(trigit)-FS*0.001)));
                else
                    v(chnnum(chnit),trigN(trigit):size(v,2))=v(chnnum(chnit),trigN(trigit):size(v,2))-(mean(v(chnnum(chnit),trigN(trigit)+FS*0.005:trigN(trigit)+FS*0.025))-mean(v(chnnum(chnit),trigN(trigit)-FS*0.020:trigN(trigit)-FS*0.001)));
                end
                v(chnnum(chnit),trigN(trigit)-FS*0.001:trigN(trigit)+FS*0.002)=interpolate(v(chnnum(chnit),trigN(trigit)-FS*0.001:trigN(trigit)+FS*0.002),1);%%%blank artefact
            end
            if N==1
                thresh=1;
            end
            [Sp,Spktimes,thresh]=OnlineSpkExtract(v(chnnum(chnit),:),thresh);
            Spk_array{chnit}=[Spk_array{chnit}; Spktimes+(N-1)*T*1000, Sp];
        end
    end
end
fclose all;
%process max four current levels - find max/min and select closest to 0.25
%and 0.75
%min 2 current levels
[Spike_trialstruct,baslinespike_trialstruct] = OnlinesortTrials(trig,Spk_array,chnnum);

%% plot sigmoid -  this can only do one channel at the moment
%rejig so record electrode just loops through the outside
loadAMP_all;
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
chn=unique(trialinfo(:,2));
chn(chn==0)=[];
count=0;
allsingletrialIDs=[];
IDs=cell(length(chn),1);

for trialloop=1:numsim_elect:length(trialinfo(:,2))
    if sum(trialinfo(trialloop:trialloop+numsim_elect-1,2)==0)==numsim_elect-1
        count=count+1;
        for chnloop=1:length(chn)
            if any(trialinfo(trialloop:trialloop+numsim_elect-1,2)==chn(chnloop))
                IDs{chnloop}(count,1:2)=[trialinfo(trialloop,17) trialloop];
            end
        end
        %allsingletrialIDs(count)=trialloop;%index   trialinfo(trialloop,1);
    end
end
        tparams = dir('*_exp_datafile_*.mat');
        tparams = tparams.name;
        load(tparams,'simultaneous_stim');
spk_all=[];
for recordchn=1:length(chnnum)
    for stimchn=1:length(chn)
        IDs{stimchn}=sortrows(IDs{stimchn});
        IDs{stimchn}=IDs{stimchn}(any(IDs{stimchn},2),:);

        for ID=1:size(IDs{stimchn},1)
            t=['ID_' num2str((IDs{stimchn}(ID,2)))];
            spkcount=zeros(40,1);
            for j=1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t),'AsArray',true),2)
                t2=['Trial_' num2str(j)];
                spkcount(j)=size(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t).(t2),1);
            end
            spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))])(ID)=mean(spkcount(1:size(struct2table(Spike_trialstruct.(['Chn_' num2str(chnnum(recordchn)-1)]).(t),'AsArray',true),2)))/6*1000;
        end
    end
end

%toc
if length(chnnum)>16
    chnpos=Depth(1);
    figure
    %%plot sigmoid
    sz=ceil(length(chnnum)/16);
    for recordchn=1:length(chnnum)
        subplot(16,sz,find(chnpos==chnnum(recordchn)))
        hold on
        for stimchn=1:length(chn)
            plot(IDs{stimchn}(:,1),spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]))
        end
        xlabel('Current (\muA)')
        ylabel('Sp/s')
    end
else
    %%plot sigmoid
    for recordchn=1:length(chnnum)
        figure
        hold on
        for stimchn=1:length(chn)
            plot(IDs{stimchn}(:,1),spk_all.(['Chn_' num2str(chnnum(recordchn)-1)]).(['stimchn_' num2str(chn(stimchn))]))
        end
        xlabel('Current (\muA)')
        ylabel('Sp/s')
        title(chnsplit{recordchn})
    end
end

%% raster
% ID=34; % there will be an increase at 100ms due to offset
% for recordchn=1:length(chnnum)
% Online_raster(ID,Spk_array{recordchn})
% end

end