penetrationDateTimeEnd=[220809211100 220811062400 220907171200 220908133500 221011133500 221012104500 221122160000 221122213000 221123111500];

%multielectrode processing
%% loop through data
D_data=dir;
namedat={D_data.name};
checkfolder=false(length(D_data),1);
IDstructsave=cell(length(D_data),1);
savefilename=cell(length(D_data),1);
parfor k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    try
        cd([D_data(k).folder filesep currD])
        filepath=pwd;
        [filepathm,name,ext] = fileparts(filepath);
        %amp(k,1:5)=loadAMP;
        sp=load([name(1:end-14) '.sp.mat'],'sp');
        sp=sp.sp;
        checkfolder(k)=true;
    catch
        continue
    end
    [stimChn,~]=loadstimulationchannels;
    chn_range=1:128;
    if str2double(name(end-12:end-7))<220812 || str2double(name(end-5:end))==072832 || str2double(name(end-5:end))==084009
        warning('port D bad')
        chn_range(chn_range>96)=[];
    end
    timestart=55;
    timeend=90;% time to stop looking for spikes in ms
    trig=loadTrig(0);
    TrialParams=loadTrialParams;
    maxid=max(cell2mat(TrialParams(:,2)));
    endtrial=maxid;
    try
        [IDstruct, baslinespikestruct,ratestruct] = sortTrials_SM(timestart,timeend,trig,0,1,1,endtrial,chn_range);
    catch
        continue
    end
    stimfilename=dir('*exp_datafile_*');
    IDstructsave{k}={IDstruct baslinespikestruct ratestruct};
    stimVar=load(stimfilename.name,'AMP','CHN');
    ti=loadTrialInfo;
    numelect=find(diff(cell2mat(ti(2:end,1)))~=0,1,'first');
    if numelect==2
        TrialsSingleElect=cell2mat(ti([false; cell2mat(ti(2:end,18))==-1],1)) ;
        ampsingle=zeros(length(TrialsSingleElect),1);
        chnsingle=zeros(length(TrialsSingleElect),1);
        pulsescountsingle=zeros(length(TrialsSingleElect),1);
        freqsingle=zeros(length(TrialsSingleElect),1);
        for i = 1:length(TrialsSingleElect)-numelect
            ampsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],18));
            chnsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],2));
            pulsescountsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],10));
            freqsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],11));
        end
    elseif numelect>2
        TrialsSingleElect=cell2mat(ti([false; cell2mat(ti(2:end,18))==-1],1)) ;
        % Find unique numbers in the input array
        unique_numbers = unique(TrialsSingleElect);
        
        % Count the occurrences of each unique number
        occurrences = histc(TrialsSingleElect, unique_numbers);
        TrialsSingleElect=unique_numbers(occurrences==(numelect-1));
        ampsingle=zeros(length(TrialsSingleElect),1);
        chnsingle=zeros(length(TrialsSingleElect),1);
         freqsingle=zeros(length(TrialsSingleElect),1);
         pulsescountsingle=zeros(length(TrialsSingleElect),1);
        for i = 1:length(TrialsSingleElect)
            ampsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],18));
            chnsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],2));
            pulsescountsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],10));
            freqsingle(i)=cell2mat(ti([false; cell2mat(ti(2:end,1))==TrialsSingleElect(i) & cell2mat(ti(2:end,18))~=-1],11));
        end
    else
        TrialsSingleElect=transpose(1:size(ti,1)-1);
        ampsingle=cell2mat(ti(2:end,18));
        chnsingle=cell2mat(ti(2:end,2));
        pulsescountsingle=cell2mat(ti(2:end,10));
        freqsingle=cell2mat(ti(2:end,11));
    end
    savefilename{k}=[{str2double(stimfilename.name(end-6:end-4))} {stimVar} {stimfilename.folder(end-12:end-7)} {[TrialsSingleElect ampsingle chnsingle pulsescountsingle freqsingle]} {stimfilename.folder(end-12:end)}];
    
end
%folders=D_data(checkfolder,:);

IDstructsave=IDstructsave(checkfolder,:);
savefilename=savefilename(checkfolder);
IDstructsave(cellfun(@isempty,IDstructsave))=[];
savefilename(cellfun(@isempty,savefilename))=[];
%% Concat trials with same stimulus file
%filenameunique=unique(savefilename);
checkoff=false(length(savefilename),1);
IDstructcompiled=cell(length(IDstructsave),1);
BaselineIDstructcompiled=cell(length(IDstructsave),1);
ratesructcompiled=cell(length(IDstructsave),1);
storefilechn=zeros(length(IDstructsave),24);
storefileamp=zeros(length(IDstructsave),2040);
storefilepulse=zeros(length(IDstructsave),500);
storefiledate=zeros(length(IDstructsave),1);
savefilenamecompiled=cell(length(IDstructsave),1);
uniqueamp=0;
for loop=1:length(IDstructsave)% loop through the stimulation pairs. Avoid using the first ones
    trialend=length(fieldnames(IDstructsave{loop}{1}));
    filessame=all(storefilechn(:,1:length(savefilename{loop}{2}.CHN(:)))==savefilename{loop}{2}.CHN(:)',2) & all(storefileamp(:,1:length(savefilename{loop}{2}.AMP(:)))==savefilename{loop}{2}.AMP(:)',2)...
        & all(storefilepulse(:,1:length(savefilename{loop}{1,4}(:,4)))==savefilename{loop}{1,4}(:,4)',2) & (storefiledate+1)>=(str2double(savefilename{loop}{3})) & (storefiledate-1)<=(str2double(savefilename{loop}{3}));
    if ~any(filessame)
        uniqueamp=uniqueamp+1;
        storefileamp(uniqueamp,1:length(savefilename{loop}{2}.AMP(:)))=savefilename{loop}{2}.AMP(:)';
        storefilechn(uniqueamp,1:length(savefilename{loop}{2}.CHN(:)))=savefilename{loop}{2}.CHN(:)';
        storefilepulse(uniqueamp,1:length(savefilename{loop}{1,4}(:,4)))=savefilename{loop}{1,4}(:,4)';
        storefiledate(uniqueamp,:)=str2double(savefilename{loop}{3});
        for trial=1:trialend
            tnum=['T', num2str(trial)];
            %if ~isempty(IDstructsave{loop}{1}.(tnum))
                numnantoadd=128-size(IDstructsave{loop}{1}.(tnum),1);
                ids_full=[IDstructsave{loop}{1}.(tnum); nan(numnantoadd,size(IDstructsave{loop}{1}.(tnum),2))];
                bss_full=[IDstructsave{loop}{2}.(tnum); nan(numnantoadd,size(IDstructsave{loop}{2}.(tnum),2))];
                rate_full=[IDstructsave{loop}{3}{trial};nan(numnantoadd,size(IDstructsave{loop}{3}{trial},2),size(IDstructsave{loop}{3}{trial},3))];
                IDstructcompiled{uniqueamp}.(tnum)=ids_full;
                BaselineIDstructcompiled{uniqueamp}.(tnum)=bss_full;
                ratesructcompiled{uniqueamp}.(tnum)=rate_full;
                savefilenamecompiled{uniqueamp}=savefilename{loop};
                
           % end
        end
    else
        index=find(filessame);
        for trial=1:trialend
            tnum=['T', num2str(trial)];
            if ~isempty(IDstructsave{loop}{1}.(tnum))
                numnantoadd=128-size(IDstructsave{loop}{1}.(tnum),1);
                ids_full=[IDstructsave{loop}{1}.(tnum); nan(numnantoadd,size(IDstructsave{loop}{1}.(tnum),2))];
                bss_full=[IDstructsave{loop}{2}.(tnum); nan(numnantoadd,size(IDstructsave{loop}{2}.(tnum),2))];
                rate_full=[IDstructsave{loop}{3}{trial};nan(numnantoadd,size(IDstructsave{loop}{3}{trial},2),size(IDstructsave{loop}{3}{trial},3))];
                IDstructcompiled{index}.(tnum)=[IDstructcompiled{index}.(tnum) ids_full];
                BaselineIDstructcompiled{index}.(tnum)=[BaselineIDstructcompiled{index}.(tnum) bss_full];
                if isempty(ratesructcompiled{index}.(tnum))
                    ratesructcompiled{index}.(tnum)=rate_full;
                else
                ratesructcompiled{index}.(tnum)=cat(3,ratesructcompiled{index}.(tnum), rate_full);
                end
            end
        end
    end
end
% remove dud channel
% trial=21:25;
% channel=105;
% file='E:\DATA\CJ_V1sigmoid\sigmoidV1E1E8_221010_205055';
% E:\DATA\CJ_V1sigmoid\PEN3_V1stim_221122_231641 - chn 108, t 25
% folder=23;
% for trials=trial
%     tnum=['T' num2str(trials)];
% IDstructcompiled{folder}.(tnum)(channel,:)=nan;
% BaselineIDstructcompiled{folder}.(tnum)(channel,:)=nan;
% ratesructcompiled{folder}.(tnum)(channel,:,:)=nan;
% end
%%
IDstructsavecompiled=[IDstructcompiled BaselineIDstructcompiled ratesructcompiled];
%%
% REMOVE BAD CHANNELS V1 (and V2 but not as affected)
numfolderstotal=size(IDstructsavecompiled,1)-sum(cellfun(@isempty, IDstructsavecompiled));
chnstoremove=[];
for loop=1:numfolderstotal(1)% loop through the stimulation pairs. Avoid using the first ones
    for trial=1:length(fieldnames(IDstructsavecompiled{loop,3}))
        tnum=['T' num2str(trial)];
        for chn=1:128
            %remove bad artefact data with baseline check
            data=squeeze(IDstructsavecompiled{loop,3}.(tnum)(chn,1:85,:));
            dat=movmean(data,3,1);
            thresh=1.5;
            checkdat=dat>thresh;
            test=sum(checkdat,2);
            numtrials=size(dat,2);
            %remove baseline data with post stim check
            dat2=squeeze(IDstructsavecompiled{loop,3}.(tnum)(chn,110:170,:));
            %dat2=movmean(data2,2,1);
            checkdat=dat2>=thresh;
            test2=sum(checkdat,2);
            checkdat2=dat2>=1;
             test3=sum(checkdat2,2);
             dat3=movmean(dat2,2,1);
            test4=sum(dat3,2);
            if any(test>numtrials/3) || any(test2>numtrials/2) || any(test3>numtrials-5) || any(test4>26)% || any(test2>5 & test3>=numtrials/2 )%remove those with bad artefacts
                %find those with bad artefacts via the baseline and then
                %remove all trials with that channel - some artefacts are
                %also bad from same channels but it happens after baseline
                %so hard to reject singularly
                chnstoremove=[chnstoremove chn];
                IDstructsavecompiled{loop,3}.(tnum)(chn,:,1:numtrials)=nan;
                IDstructsavecompiled{loop,2}.(tnum)(chn,1:numtrials)=nan;
                IDstructsavecompiled{loop,1}.(tnum)(chn,1:numtrials)=nan;
                %                 for ttemp=1:length(fieldnames(IDstructsavecompiled{loop,3}))
                %
                % %                     tnumtemp=['T' num2str(ttemp)];
                % %                     IDstructsavecompiled{loop,3}.(tnumtemp)(chn,:,1:numtrials)=nan;
                % %                     IDstructsavecompiled{loop,2}.(tnumtemp)(chn,1:numtrials)=nan;
                % %                     IDstructsavecompiled{loop,1}.(tnumtemp)(chn,1:numtrials)=nan;
                %                 end
            end
            dat2=dat>0;
            periodic50hz=20;%1 every 20ms
            for i=1:periodic50hz %remove periodic noise in baseline trials
                check=i:periodic50hz:85;
                check2=i+periodic50hz/2:periodic50hz:85;
                removetrials=all(dat2(check,:)==1) & all(dat2(check2,:)==0);
                if sum(removetrials==1)>2%its more likely to have multiple trials if periodic
                    IDstructsavecompiled{loop,3}.(tnum)(chn,:,1:numtrials)=nan;
                    IDstructsavecompiled{loop,2}.(tnum)(chn,1:numtrials)=nan;
                    IDstructsavecompiled{loop,1}.(tnum)(chn,1:numtrials)=nan;
                end
            end
        end
    end
end
%% remove bursty channels
timebefore=0.09;%in seconds 0.09
timeafter=0.09;%in seconds0.09
FS=30000;
numfolderstotal=size(IDstructsavecompiled,1)-sum(cellfun(@isempty, IDstructsavecompiled));
for loop=1:numfolderstotal(1)% loop through the stimulation pairs. Avoid using the first ones
    fold_int=dir(['*' savefilenamecompiled{loop}{5}]);
    cd(fold_int.name)
    spall=load([fold_int.name(1:end-14) '.sp.mat'],'sp');
    trig=loadTrig;
    TrialParams=loadTrialParams;
    numelect=find(diff(cell2mat(TrialParams(:,1)))~=0,1,'first');
    for chn=1:128
        sp=spall.sp{chn};
        wsp=removeArtifactSpikes(sp);
        for trial=1:length(fieldnames(IDstructsavecompiled{loop,3}))
            tnum=['T' num2str(trial)];
            TrialParamstID = find(cell2mat(TrialParams(1:numelect:end,2)) == trial); %identifies trial row matching trial ID
            trigtID = trig(TrialParamstID);
            trigtID(trigtID==-500)=[];
            
            data=squeeze(IDstructsavecompiled{loop,3}.(tnum)(chn,1:180,:));
            removetrial=[];
            
            
            if any(data>2,'all')
                dat=movmean(data,3,1);
                thresh=mean(dat,2,'omitnan')+std(dat,[],2,'omitnan').*3;
                checkdat=dat>thresh;
                
                for trialsnum=1:size(checkdat,2)
                    f = find(diff([0,checkdat(:,trialsnum)',0]==1));
                    p = f(1:2:end-1);  % Start indices
                    numconsec = f(2:2:end)-p;  % Consecutive ones� counts
                    if any(numconsec>15) || any(numconsec>10 & p>90+10 & chn>64)
                        removetrial=[removetrial trialsnum];
                    end
                end
                datsummed=sum(dat,2,'omitnan');
                ds_raw=sum(data,2,'omitnan');
                if (sum(datsummed(1:90))<5 && sum(datsummed(91:180))>500 && all(datsummed(91:180)<30))
                    %if very bad remove that channel from all folders will likely cause problems throughout dataset -
                    %should be 109 and 108
                    for loopfold=1:numfolderstotal(1)
                        if savefilenamecompiled{loopfold}{3}(1:4)==savefilenamecompiled{loop}{3}(1:4)
                            for ittrial=1:length(fieldnames(IDstructsavecompiled{loopfold,3}))
                                trialnum=['T' num2str(ittrial)];
                                IDstructsavecompiled{loopfold,3}.(trialnum)(chn,:,:)=nan;
                                IDstructsavecompiled{loopfold,2}.(trialnum)(chn,:)=nan;
                                IDstructsavecompiled{loopfold,1}.(trialnum)(chn,:)=nan;
                            end
                        end
                    end
                    break;
                end
                
                
                %remove channels with bad artifact
                number_artifacts=0;
                number_artifacts_2=0;
                for indT=1:length(trigtID)
                    offsetspk=trigtID(indT)*1000./FS;
                    spktimes=sp(((sp(:,1)>(offsetspk-timebefore*1000))&(sp(:,1)<(offsetspk+timeafter*1000))),:);
                    number_artifacts=number_artifacts+sum(spktimes(:,2:end)<-200,'all')+sum(spktimes(:,2:end)>200,'all');
                    spktimes_2=sp(((sp(:,1)>(offsetspk))&(sp(:,1)<(offsetspk+0.02*1000))),:);
                    number_artifacts_2=number_artifacts_2+sum(spktimes_2(:,2:6)>50,'all');
                    if (number_artifacts>5 || number_artifacts_2>5) && (any(datsummed(112:118)>15)|| any(datsummed(150:155)>15) || any(ds_raw(90:110)>9 & sum(ds_raw(1:89))<6) || any(ds_raw(90:110)>sum(ds_raw(1:89))*3)||any(ds_raw(140:150)>sum(ds_raw(1:89))))
                        removetrial=1:size(checkdat,2);
                        %                             figure
                        %                             plot(1*1000/FS:1000/FS:49*1000/FS,spktimes(:,2:end)')
                        %                             xlabel('Time (ms)')
                        %                             ylabel('Spike amplitude (uV)')
                        %                             title(['fold ' num2str(loop) ' chn ' num2str(chn) ' trial ' num2str(trial)])
                        break
                    end
                    spktimes=wsp(((wsp>(offsetspk-timebefore*1000))&(wsp<(offsetspk+timeafter*1000))),:)-offsetspk;
                    if ~isempty(spktimes)
                        rate_wsp=hist(spktimes,-timebefore*1000:timeafter*1000);
                        IDstructsavecompiled{loop,3}.(tnum)(chn,1:timeafter*1000+timebefore*1000+1,indT)=squeeze(IDstructsavecompiled{loop,3}.(tnum)(chn,1:timeafter*1000+timebefore*1000+1,indT))-rate_wsp;
                    end
                end
            end
            
            
            if ~isempty(removetrial)
                IDstructsavecompiled{loop,3}.(tnum)(chn,:,removetrial)=nan;
                IDstructsavecompiled{loop,2}.(tnum)(chn,removetrial)=nan;
                IDstructsavecompiled{loop,1}.(tnum)(chn,removetrial)=nan;
            end
          
        end
    end
    cd(fold_int.folder)
end
%% for multipulse
%remove chn 108 in  C:\data\multipulse\PEN3_V1stim_2pulse_221123_074018
for trial=1:81
    tnum=['T' num2str(trial)];
    IDstructsavecompiled{7,1}.(tnum)(108,:,:)=nan;
    IDstructsavecompiled{8,1}.(tnum)(111,:,:)=nan;
end
multipulseplotdata(IDstructsavecompiled(:,3),savefilenamecompiled,'v1')
%% plot spiking rate epoch

excitesupress=1;%1 for excite, 0 for supress - lower threshold for supress(see below)
chnrange=1:64;
normalisedat=0;
close all
if all(savefilenamecompiled{7,1}{1,2}.CHN<65)
    chnstoremove=[];
end
epochratespike(IDstructsavecompiled(:,3),savefilenamecompiled,chnrange,excitesupress,chnstoremove,normalisedat);
%% multielectrode plot
multielectplot(IDstructsavecompiled(:,3),savefilenamecompiled,chnrange,excitesupress,chnstoremove)
%% Overlap RFs
chnrange=1:64;
%change chnrange for V1 vs V2 then use these below to select overlap amount
fulloverlapRFs=false;%only those with full overlap
NOTfulloverlapRFs=true;%everything without full overlap

%don't touch
partialoverlapRFs=false; %only partial overlap
AlloverlapRFs=false;%find stim data from overlapping RFs = true or false for non-overlap

excitesupress=1;%1 for excite, 0 for supress - lower threshold for supress(see below)

normalisedat=0;
close all
% if all(savefilenamecompiled{3,1}{1,2}.CHN<65)
%     chnstoremove=[];
% end

%CJ219 had no overlap, CJ222 had partial overlap both pens, CJ225 full 
%overlap both pens, CJ230 had partial overlap in p3 and full overlap p1
date_time_partialoverlap=[{'220906_230000'} {'220907_180000'};{'220907_190000'} {'220908_140000'}; {'221122_220000'} {'221123_120000'};];%start and end penetration folders
date_time_fulloverlap=[{'221010_180000'} {'221012_110000'};{'221121_180000'} {'221122_160000'}; ];%start and end penetration folders
if fulloverlapRFs || NOTfulloverlapRFs
    date_time=date_time_fulloverlap;
elseif partialoverlapRFs
    date_time=date_time_partialoverlap;
else
    date_time=[date_time_partialoverlap;date_time_fulloverlap];
end

%get stim data from relevent penetration
savefilenamecompiled=savefilenamecompiled(~cellfun(@isempty, savefilenamecompiled));
IDstructsavecompiled=IDstructsavecompiled(~cellfun(@isempty, IDstructsavecompiled(:,1)),1:3);
datetimeStamps=cell(1,length(savefilenamecompiled));
for foldit=1:length(savefilenamecompiled)
        datetimeStamps{foldit}=savefilenamecompiled{foldit}{5};
end
% Convert datetime stamps to datetime objects
formats = 'yyMMdd_HHmmss';  % Format of your datetime stamps
datetimes = datetime(datetimeStamps, 'InputFormat', formats);

invertRFrange=false;
withinrangeAllpens=false(1,length(datetimes));
for pen_times=1:size(date_time,1)
% Define your datetime range
startRange = datetime(date_time{pen_times,1}, 'InputFormat', formats);
endRange = datetime(date_time{pen_times,2}, 'InputFormat', formats);

% Check if each datetime falls within the range
 isWithinRange = datetimes >= startRange & datetimes <= endRange;
 
if NOTfulloverlapRFs || (~AlloverlapRFs && ~fulloverlapRFs && ~partialoverlapRFs)
   invertRFrange=true;
end
withinrangeAllpens=withinrangeAllpens|isWithinRange;
end
if invertRFrange
    withinrangeAllpens=~withinrangeAllpens;
end
% get relevent folders from penetration
relevantData=IDstructsavecompiled(withinrangeAllpens,1:3);
relevantSavefilename=savefilenamecompiled(withinrangeAllpens);

alldata=epochratespike(relevantData(:,3),relevantSavefilename,chnrange,excitesupress,chnstoremove,normalisedat);
%% multielectrode plot
multielectplot(relevantData(:,3),relevantSavefilename,chnrange,excitesupress,chnstoremove)
%%
fulloverlapdata=alldata;
noneApartial=alldata;
%% overlap,non-overlapping
V2=true;%false if V1
NPO_data=[noneApartial fulloverlapdata];%full overlap or no overlap

s = RandStream('mt19937ar','Seed',296);
subselectbaseline=1:89;
permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
subselectbaseline=permuted_numbers(1:50);

if V2
    timetoavg=92:181;
else
    timetoavg=92:112;
end

baselinenpo=cellfun(@(x) mean(mean(x(:,1:89)-mean(x(:,subselectbaseline),2,'omitnan'),2,'omitnan')),NPO_data,'UniformOutput',true).*1000;
baselinestdnpo=cellfun(@(x) mean(SEM(x(:,1:89)-mean(x(:,subselectbaseline),2,'omitnan'),0)),NPO_data,'UniformOutput',true).*1000;
avgnpo=cellfun(@(x) mean(mean(x(:,timetoavg)-mean(x(:,subselectbaseline),2,'omitnan'),2,'omitnan')),NPO_data,'UniformOutput',true).*1000;

stdnpo=cellfun(@(x) mean(SEM(x(:,timetoavg)-mean(x(:,subselectbaseline),2,'omitnan'),0)),NPO_data,'UniformOutput',true).*1000;
figure;
errorbar([0 2 5 6 8 10],[(baselinenpo(5,1)); avgnpo(:,1)],[(baselinestdnpo(5,1)); stdnpo(:,1)],'r')
hold on
errorbar([0 2 5 6 8 10],[(baselinenpo(5,2)); avgnpo(:,2)],[(baselinestdnpo(5,2)); stdnpo(:,2)],'b')
%errorbar([0 2 5 6 8 10],[mean(baselinenpo(:,3)); avgnpo(:,3)],[mean(baselinestdnpo(:,3)); stdnpo(:,3)],'g')
set(gca,'TickDir','out');
xlabel('Current (\muA)');
ylabel('Average Firing Rate (Sp/s)')
if V2
title('V2')
else
    title('V1')
end
legend('Out RF','In RF')
beautifyPlot;

figure; 
ax=axes;
hold on;
stdshade((NPO_data{5,1}).*1000,0.2,'r',[-90:90],1,ax);
stdshade((NPO_data{5,2}).*1000,0.2,'b',[-90:90],1,ax);
xlabel('Time(ms)')
ylabel('FR (sp/s)')
if V2
title('V2')
else
    title('V1')
end
 xlim([-50,89])
set(gca,'TickDir','out');
beautifyPlot;

% big small
avgallstim=cellfun(@(x) sort((mean(x(:,timetoavg)-mean(x(:,subselectbaseline),2,'omitnan'),2,'omitnan')).*1000),NPO_data,'UniformOutput',false);
smalldat=cellfun(@(x) mean(x(1:round(length(x)/2))),avgallstim,'UniformOutput',true);
bigdat=cellfun(@(x) mean(x(round(length(x)/2):end)),avgallstim,'UniformOutput',true);

figure
plot(smalldat)
hold on
plot(bigdat)


%%
%manual input required for penetration info
chnrange=1:64;%pick V1 or V2 channels
overlappingRFs=1;%1 for overlapping channels only, -1 for non overlapping, 0 for all channels
if overlappingRFs~=0
    load('E:\DATA\MarmoblitzVisualDataCJ\2022\09\07\RF_MAP_full.mat')
end
date_time=[{'220906_220000'} {'220907_180000'}];%start and end penetration folders

%automatic from here

%find overlapping RFs - use threshold 95% of peak? 
% Step 1: Find the peak location in each 20x20 double
peak_locations = cell(1, 128);
for chn = 1:128
    data = RF{chn};
    [max_value, max_index] = max(data(:));
    [row, col] = ind2sub(size(data), max_index);
    peak_locations{chn} = [row, col];
end

% Step 2: Determine the most common peak location
all_peak_locations = cell2mat(peak_locations');
location_counts = accumarray(all_peak_locations, 1);
[rf_R,rf_C]=find(mean(location_counts,'all')+std(location_counts,[],'all')*3<location_counts);%has to be 3sd above mean to be the RF

% Step 3: Split the cell array into two groups based on the peak location
inRF = {};
outRF = {};
chn_inRF=[];
chn_outRF=[];

for chn = 1:128
    if any(peak_locations{chn}(1)==rf_R & peak_locations{chn}(2)==rf_C)
        inRF{end+1} = RF{chn};%in RF
        chn_inRF=[chn_inRF chn];
    else
        outRF{end+1} = RF{chn};%not in RF
        chn_outRF=[chn_outRF chn];
    end
end

%get stim data from relevent penetration
savefilenamecompiled=savefilenamecompiled(~cellfun(@isempty, savefilenamecompiled));
IDstructsavecompiled=IDstructsavecompiled(~cellfun(@isempty, IDstructsavecompiled(:,1)),1:3);
datetimeStamps=cell(1,length(savefilenamecompiled));
for foldit=1:length(savefilenamecompiled)
        datetimeStamps{foldit}=savefilenamecompiled{foldit}{5};
end
% Convert datetime stamps to datetime objects
formats = 'yyMMdd_HHmmss';  % Format of your datetime stamps
datetimes = datetime(datetimeStamps, 'InputFormat', formats);

% Define your datetime range
startRange = datetime(date_time{1}, 'InputFormat', formats);
endRange = datetime(date_time{2}, 'InputFormat', formats);

% Check if each datetime falls within the range
isWithinRange = datetimes >= startRange & datetimes <= endRange;

% get relevent folders from penetration
relevantData=IDstructsavecompiled(isWithinRange,1:3);
relevantSavefilename=savefilenamecompiled(isWithinRange);

% now need to run epoch spike rate with relevant V1 and V2 channels found
% in chn_inRF VS chn_outRF
if overlappingRFs==1
    chns=chn_inRF;
elseif overlappingRFs==-1
    chns=chn_outRF;
end


%CJ225
% V1
% S1E 2 3 7 8 
% S4E whole shank, s11-16 weak 8-10 bursty
% S2E 2 3 4 5 6 7 8 9 10
% S3E 1-9
% 
% V2
% S1E 3 4 5 6 7 8 9 10 11 13 15 16
% S4E 2 3 8 9 10
% S2E 3 4 6 10 12 13
% S3E 5 7 8 9 10

chs_inRFCJ225=[1 4 90 87 35 64 113 100];
chn_outRFCJ225=setdiff(1:128, chs_inRFCJ225);

date_time=[{'221011_160000'} {'221012_120000'}];%start and end penetration folders
% Define your datetime range
startRange = datetime(date_time{1}, 'InputFormat', formats);
endRange = datetime(date_time{2}, 'InputFormat', formats);

% Check if each datetime falls within the range
isWithinRange = datetimes >= startRange & datetimes <= endRange;
relevantDataCJ225=IDstructsavecompiled(isWithinRange,1:3);
relevantSavefilenameCJ225=savefilenamecompiled(isWithinRange);

relevantData_all=[relevantData;relevantDataCJ225];
savefilename_all=[relevantSavefilename;relevantSavefilenameCJ225];

if overlappingRFs==1
    chnsCJ225=chs_inRFCJ225;
elseif overlappingRFs==-1
    chnsCJ225=chn_outRFCJ225;
end

%make the channel lists the same size to ensure we can concatenate them
largestRFsize=max(length(chns),length(chnsCJ225));
chnsCJ225=[chnsCJ225 zeros(1,largestRFsize-length(chnsCJ225))];
chns=[chns zeros(1,largestRFsize-length(chns))];

chnCJ222=repmat(chns,[length(relevantData(:,3)),1]);
chnCJ225=repmat(chnsCJ225,[length(relevantDataCJ225(:,3)),1]);

chns_all=[chnCJ222;chnCJ225];

excitesupress=1;%1 for excite, 0 for supress - lower threshold for supress(see below)

normalisedat=0;
close all

%CJ230 - overlapping in pen1
%overlapping in pen2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%messed with significance - need to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%check this
epochratespikeRFmatch(relevantData_all(:,3),savefilename_all,chnrange,excitesupress,chnstoremove,normalisedat,chns_all);






%epochratespikeRFmatch(relevantData(:,3),relevantSavefilename,chnrange,excitesupress,chnstoremove,normalisedat,chns);

%%
                    
E_MAP=Depth(1);
[mappedchns,~]=find(E_MAP==chns);
intersect(relevantSavefilename{4}{1, 2}.CHN,mappedchns)
intersect([67 86 116 99],mappedchns)


%%
close all
excitesupress=1;%1 for excite, 0 for supress - lower threshold for supress(see below)
chnrange=1:64;
multielectplot(IDstructsavecompiled(:,3),savefilenamecompiled,chnrange,excitesupress,chnstoremove)

%%
relationshipV1V2chn(IDstructsavecompiled(:,3),savefilenamecompiled,chnstoremove) %%%%%%%continue working on this function

%% V2 based on V1 resp

%% plot rate from individual penetration
%E:\DATA\CJ_V1sigmoid\Pen1_E6E13_221121_220434 - early resp - folder 13
folder=11;
trial=5;
plotV1V2matrix(IDstructsavecompiled{folder,3},trial)



%% trial average
storefileamp(all(storefileamp==0,2),:)=[];
storefilechn(all(storefilechn==0,2),:)=[];

cd([D_data(4).folder filesep D_data(4).name;])
avgspiking=cell(size(storefileamp,1),1);
ratespiking=cell(size(storefileamp,1),1);
for numerfolders=1:size(storefileamp,1)
    [avgnospT,stderrspktrial,~]=AverageTrialResponse_SM(IDstructcompiled{numerfolders},BaselineIDstructcompiled{numerfolders});
    avgspiking{numerfolders}=avgnospT;
    for trial=1:trialend
        tnum=['T', num2str(trial)];
        ratespiking{numerfolders}(:,:,trial)=nanmean(ratesructcompiled{numerfolders}.(tnum),3);
    end
end

%% sort into responding not responding using standard deviations
ampinterest=6;
ratespiking_nill=cell(size(storefileamp,1),1);
ratespiking_excite=cell(size(storefileamp,1),1);
ratespiking_supress=cell(size(storefileamp,1),1);
for numfolders=1:size(storefileamp,1)
    ampit=reshape(storefileamp(numfolders,:),[],size(storefilechn,2));%%%%%%check this!!! need to see if moving has affected
    if any(ampit~=-1 & ampit~=ampinterest,'all')% pull out 6uA results
        continue
    end
    for chn=1:64
        tmp=squeeze(ratespiking{numfolders}(chn,:,255));
        %         epoch=10;%ms
        %         TimeWindow=0;
        %         for time=0:floor(90/epoch)-2
        %             if  ttest(tmp(75:85),tmp(92+(time*epoch):92+((time+1)*epoch)))==1
        %                 TimeWindow=92+(time*epoch):92+((time+1)*epoch);
        %                 break
        %             end
        %         end
        
        if any((mean(tmp(1:85))+std(tmp(1:85))*3)<tmp(92+0:92+8))%92:92+84
            ratespiking_excite{numfolders}=cat(1,ratespiking_excite{numfolders}, ratespiking{numfolders}(chn,:,:));
        elseif any((mean(tmp(1:85))-std(tmp(1:85)))>tmp(92:92+84))
            ratespiking_supress{numfolders}=cat(1,ratespiking_supress{numfolders}, ratespiking{numfolders}(chn,:,:));
        else
            ratespiking_nill{numfolders}=cat(1,ratespiking_nill{numfolders}, ratespiking{numfolders}(chn,:,:));
        end
    end
end




%% trial avg single elect - - standard deviation - only one peak 09/05/23 - bins all electrodes based on thresh
excitesupress=1;%1 for excite, 0 for supress - lower threshold for supress(see below)
cd([D_data(6).folder filesep D_data(6).name;])
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,[1,3,4,2]);
alldata=cell(length(savefilename{1}{2}.AMP),1);
trialsig=cell(9,length(savefilename{1}{2}.AMP),length(savefilename));
folderIDsig=cell(9,length(savefilename{1}{2}.AMP));
stimchncount=0;
spread=nan(length(savefilename{1}{2}.AMP),500);
spreadgroup=nan(length(savefilename{1}{2}.AMP),9,500);
groupdata=cell(length(savefilename{1}{2}.AMP),1);
for ampit=1:length(savefilename{1}{2}.AMP)
ampinterest=savefilename{1}{2}.AMP(ampit);
groupdata{ampit}=cell(9,1);
significantspread=nan(128,500);
totalchncount=0;
significantspread_groupsplit=nan(128,500,9);
for numerfolders=1:size(IDstructsave,1)
    trialend=length(fieldnames(IDstructsave{numerfolders}{1}));
    for trial=1:trialend
        ratespiking{numerfolders}(:,:,trial)=nanmean(IDstructsave{numerfolders}{3}{trial},3);
    end
   %%%%%%need to see if this line will pick the correct spiking trials.
   %%%%%%compare amp and select trial with that amp from savefile name.
   %%%%%%need to see if this owrrks with both the current steered data and
   %%%%%%non-current steering
    rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1));

    totalchncount=64*size(rateAMPua,3)+totalchncount;
    meancatdata=cat(2,mean(rateAMPua(:,92:101,:),2),mean(rateAMPua(:,102:111,:),2),mean(rateAMPua(:,112:121,:),2),mean(rateAMPua(:,122:131,:),2),mean(rateAMPua(:,132:141,:),2),mean(rateAMPua(:,142:151,:),2),mean(rateAMPua(:,152:161,:),2),mean(rateAMPua(:,162:171,:),2),mean(rateAMPua(:,172:181,:),2));
    sig=false(128,size(meancatdata,3),9);
    if excitesupress==0
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)-(std(rateAMPua(:,1:85,:),[],2).*1)); %%%%%NOTE LOWER THRESHOLD
        for groupit=1:9
            sig(:,:,groupit)=squeeze(meancatdata(:,groupit,:))<thresh;
        end
        check_multiple=sum(sig,3)>1;
        [~,pmax]=min(meancatdata,[],2);%pmin
        pmax=squeeze(pmax);%min
    else
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)+(std(rateAMPua(:,1:85,:),[],2).*3));
        for groupit=1:9
            sig(:,:,groupit)=squeeze(meancatdata(:,groupit,:))>thresh;
        end
        check_multiple=sum(sig,3)>1;
        [~,pmax]=max(meancatdata,[],2);
        pmax=squeeze(pmax);
    end
    

  

    for chn=1:64
        for stimchn=1:size(sig,2)
            if check_multiple(chn,stimchn)
            sig(chn,stimchn,:)=false;
             sig(chn,stimchn,pmax(chn,stimchn))=true;
            end
            group=find(squeeze(sig(chn,stimchn,:)==1));
            if ~isempty(group)
                groupdata{ampit}{group}=[groupdata{ampit}{group}; rateAMPua(chn,:,stimchn)];
                alldata{ampit}=[alldata{ampit}; rateAMPua(chn,:,stimchn)];
                significantspread(chn,stimchn+stimchncount)=meancatdata(chn,group,stimchn);
                significantspread_groupsplit(chn,stimchn+stimchncount,group)=meancatdata(chn,group,stimchn);
                trialsig{group,ampit,numerfolders}=[trialsig{group,ampit,numerfolders} stimchn];%determine how many different stim electrodes creating significance
                folderIDsig{group,ampit}=[folderIDsig{group,ampit} str2double(savefilename{numerfolders}{3})];%determine how many different monkeys creating significance
            end
            
        end
    end
    stimchncount=stimchncount+size(sig,2);
    
end
stimchncount=0;
%current significance
E_MAP=Depth;
orderedsigchns=significantspread(E_MAP,:);
avgspread=nan(size(orderedsigchns,2),1);

for itstimchn=1:size(orderedsigchns,2)
    tmp=reshape(orderedsigchns(1:64,itstimchn),16,4);
    tmp=tmp(:,[1 3 4 2]);
    tmp(isnan(tmp))=0;
    xcent=sum(tmp.*[1 2 3 4],'all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
    ycent=sum(tmp.*(1:16)','all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
 
    xval=~isnan(tmp).*[1 2 3 4];
    xval(xval==0)=nan;
    yval=~isnan(tmp).*(1:16)';
    yval(yval==0)=nan;
    avgspread(itstimchn)=mean(sqrt((((yval-ycent).*50).^2)+(((xval-xcent).*200).^2)),'all','omitnan');
    
end
%group significance
orderedgroupsig=significantspread_groupsplit(E_MAP,:,:);
avgspreadwg=nan(9,500);
for group=1:9
    for itstimchn=1:size(orderedsigchns,2)
        tmp=reshape(squeeze(orderedgroupsig(1:64,itstimchn,group)),16,4);
        tmp=tmp(:,[1 3 4 2]);
        tmp(isnan(tmp))=0;
        xcent=sum(tmp.*[1 2 3 4],'all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
        ycent=sum(tmp.*(1:16)','all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
        
        xval=~isnan(tmp).*[1 2 3 4];
        xval(xval==0)=nan;
        yval=~isnan(tmp).*(1:16)';
        yval(yval==0)=nan;
        avgspreadwg(group,itstimchn)=mean(sqrt((((yval-ycent).*50).^2)+(((xval-xcent).*200).^2)),'all','omitnan');
    end
end

spreadgroup(ampit,1:9,1:length(avgspread))=avgspreadwg;
spread(ampit,1:length(avgspread))=avgspread;
figure(1000*ampit); hold on;cellfun(@(x) plot(-90:90,mean(x,1).*1000), groupdata{ampit})
color1 = linspace(0,1,groupit);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
indvsig=cellfun(@(x) size(x,1), groupdata{ampit});
text(-80,16,'# sig:')
text(-80,10,num2str(indvsig))
totalsig=sum(indvsig);
text(-80,20,['Total sig: ' num2str(totalsig) ' / ' num2str(totalchncount)])
xlabel('Time (ms)')
ylabel('Firing rate (Sp/s)')
xlim([-85 85])
ylim([0 25])
set(gca,'TickDir','out');
title([num2str(ampinterest) '\muA'])
end
%% # Sig monkeys and stim elect per time
sigmonkeys=cellfun(@(x) unique(x), folderIDsig, 'UniformOutput', false);
sigmonkeys=cellfun(@(x) sum(diff(x)>3)+1,sigmonkeys);


sigstimelect=sum((cellfun(@(x) length(unique(x)),trialsig,'UniformOutput',true)),3);


figure;ax=axes; [lineOut, fillOut] = stdshade(spread',0.2,'r',[2 5 6 8 10],1,ax);
set(gca,'TickDir','out');
ylabel('Distance from centroid (um)')
xlabel('Current (uA)')


figure;ax=axes; hold on;
for ampit=1:5
[lineOut, fillOut] = stdshade(squeeze(spreadgroup(ampit,:,:))',0.2,[0 0 ampit/5],[1:9],1,ax);
end
set(gca,'TickDir','out');
ylabel('Distance from centroid (um)')
xlabel('Group time')


figure; heatmap(sum(heatmap_centroid{1}>0,3)); set(gca,'ColorScaling','log')
figure; heatmap(sum(heatmap_centroid{5}>0,3));set(gca,'ColorScaling','log')
%% sigmoid using groupdata
sigmoidampepoch=zeros(5,9);
for ampit=1:5
    for epoch=1:9
        sigmoidampepoch(ampit,epoch)=mean(alldata10uA{ampit}{epoch}(:,82+(10*epoch):91+(10*epoch)),'all'); 
    end
end
figure; hold on
for i=1:9
plot([2 5 6 8 10],sigmoidampepoch(:,i).*1000)
end
color1 = linspace(0,1,9);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
plot([2 5 6 8 10],mean(sigmoidampepoch,2).*1000,'k')
xlabel('current (\muA)')
ylabel('Firing rate (sp/s)')
leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','61:72','72:81','82:91','Average');
title(leg,'Time epoch(ms)')
set(gca,'TickDir','out');

indvsig=zeros(5,9);
for i=1:5
indvsig(i,:)=cellfun(@(x) size(x,1), groupdata{i});
end

figure
hold on
plot([2 5 6 8 10],indvsig)
colororder(newcolors);
plot([2 5 6 8 10],sum(indvsig,2),'k')
yline((0.003)*totalchncount(1),'r')
xlabel('current (\muA)')
ylabel('# elect sig')
leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','61:72','72:81','82:91','Total','Significance');
title(leg,'Time epoch(ms)')
set(gca,'TickDir','out');




figure; hold on
for ampit=1:5
plot(mean(alldata{ampit}))
end
figure; hold on
for ampit=1:5
plot(mean(alldata{ampit}))
end

cellfun(@(x) mean(x(:,92:181),'all'),alldata10uA)

figure;
plot(mean(sigmoidampepoch,2).*1000)
%% epoch has to be increasing with current once significant
% trial avg single elect - - standard deviation -  bins all electrodes based on thresh
excitesupress=1;%1 for excite, 0 for supress - lower threshold for supress(see below)
cd([D_data(15).folder filesep D_data(15).name;])
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,[1,3,4,2]);
alldata=cell(length(savefilename{15}{2}.AMP),1);
trialsig=cell(9,length(savefilename{15}{2}.AMP),length(savefilename));
folderIDsig=cell(9,length(savefilename{15}{2}.AMP));
stimchncount=0;
spread=nan(length(savefilename{15}{2}.AMP),500);
spreadgroup=nan(length(savefilename{15}{2}.AMP),9,500);
groupdata=cell(length(savefilename{15}{2}.AMP),1);
alldata10uA=cell(5,1);
 totalchncount=zeros(5,1);
heatmap_centroid=cell(5,1);
for ampit=1:length(savefilename{15}{2}.AMP)
    heatmap_centroid{ampit}=nan(31,7,500*9);
ampinterest=savefilename{15}{2}.AMP(ampit);
groupdata{ampit}=cell(9,1);
alldata10uA{ampit}=cell(9,1);
significantspread=nan(128,500);

significantspread_groupsplit=nan(128,500,9);
for numerfolders=1:size(IDstructsave,1)
    trialend=length(fieldnames(IDstructsave{numerfolders}{1}));
    for trial=1:trialend
        ratespiking{numerfolders}(:,:,trial)=nanmean(IDstructsave{numerfolders}{3}{trial},3);
    end
   %%%%%%need to see if this line will pick the correct spiking trials.
   %%%%%%compare amp and select trial with that amp from savefile name.
   %%%%%%need to see if this owrrks with both the current steered data and
   %%%%%%non-current steering
    [m,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest)); ampinterest=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp
   if excitesupress==1
       %determines if spiking increases with increases in current overall -
       %just has to have a net increase between 2 and 10
       rateepochcurrent=nan(128,9,length(savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1)),5);
       for ampit1=1:length(savefilename{15}{2}.AMP)
           ampinterest1=savefilename{15}{2}.AMP(ampit1);
           [m,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest1)); ampinterest1=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp
           
           rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest1),1));
           rateepochcurrent(1:size(rateAMPua,1),:,:,ampit1)=cat(2,mean(rateAMPua(:,92:101,:),2),mean(rateAMPua(:,102:111,:),2),mean(rateAMPua(:,112:121,:),2),mean(rateAMPua(:,122:131,:),2),mean(rateAMPua(:,132:141,:),2),mean(rateAMPua(:,142:151,:),2),mean(rateAMPua(:,152:161,:),2),mean(rateAMPua(:,162:171,:),2),mean(rateAMPua(:,172:181,:),2));
       end
       validchannels=sum(diff(rateepochcurrent,[],4),4,'omitnan')>0;
   end
   rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1));
   totalchncount(ampit)=64*size(rateAMPua,3)+totalchncount(ampit);
   meancatdata=cat(2,mean(rateAMPua(:,92:101,:),2),mean(rateAMPua(:,102:111,:),2),mean(rateAMPua(:,112:121,:),2),mean(rateAMPua(:,122:131,:),2),mean(rateAMPua(:,132:141,:),2),mean(rateAMPua(:,142:151,:),2),mean(rateAMPua(:,152:161,:),2),mean(rateAMPua(:,162:171,:),2),mean(rateAMPua(:,172:181,:),2));
    sig=false(128,9,size(meancatdata,3));
    if excitesupress==0
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)-(std(rateAMPua(:,1:85,:),[],2).*1)); %%%%%NOTE LOWER THRESHOLD
        for groupit=1:9
            sig(:,groupit,:)=squeeze(meancatdata(:,groupit,:))<thresh;
        end
        check_multiple=sum(sig,2)>1;
        [~,pmax]=min(meancatdata,[],2);%pmin
        pmax=squeeze(pmax);%min
    else
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)+(std(rateAMPua(:,1:85,:),[],2).*3));
        for groupit=1:9
            sig(1:size(meancatdata,1),groupit,:)=squeeze(meancatdata(:,groupit,:))>thresh;
        end
        sig=sig&validchannels;
        check_multiple=squeeze(sum(sig,2)>1);
        [~,pmax]=max(meancatdata,[],2);
        pmax=squeeze(pmax);
        
    end
   

  

    for chn=1:64
        for stimchn=1:size(sig,3)
            if check_multiple(chn,stimchn)
            sig(chn,:,stimchn)=false;
             sig(chn,pmax(chn,stimchn),stimchn)=true;
            end
            group=find(squeeze(sig(chn,:,stimchn)==1));
            if ~isempty(group)
                groupdata{ampit}{group}=[groupdata{ampit}{group}; rateAMPua(chn,:,stimchn)];
                alldata{ampit}=[alldata{ampit}; rateAMPua(chn,:,stimchn)];
                if ampit==5
                for itAua=1:5
                    ampinterest1=savefilename{15}{2}.AMP(itAua);
                    [m,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest1)); ampinterest1=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp

                   rateAMPua1=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest1),1));
                    alldata10uA{itAua}{group}=[alldata10uA{itAua}{group}; rateAMPua1(chn,:,stimchn)];
                end
                end
                significantspread(chn,stimchn+stimchncount)=meancatdata(chn,group,stimchn);
                significantspread_groupsplit(chn,stimchn+stimchncount,group)=meancatdata(chn,group,stimchn);
                trialsig{group,ampit,numerfolders}=[trialsig{group,ampit,numerfolders} stimchn];%determine how many different stim electrodes creating significance
                folderIDsig{group,ampit}=[folderIDsig{group,ampit} str2double(savefilename{numerfolders}{3})];%determine how many different monkeys creating significance
            end
            
        end
    end
    stimchncount=stimchncount+size(sig,3);
    
end
stimchncount=0;
%current significance
E_MAP=Depth;
orderedsigchns=significantspread(E_MAP,:);
avgspread=nan(size(orderedsigchns,2),1);

for itstimchn=1:size(orderedsigchns,2)
    tmp=reshape(orderedsigchns(1:64,itstimchn),16,4);
    tmp=tmp(:,[1 3 4 2]);
    tmp(isnan(tmp))=0;
    xcent=sum(tmp.*[1 2 3 4],'all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
    ycent=sum(tmp.*(1:16)','all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
 
    xval=~isnan(tmp).*[1 2 3 4];
    xval(xval==0)=nan;
    yval=~isnan(tmp).*(1:16)';
    yval(yval==0)=nan;
    avgspread(itstimchn)=mean(sqrt((((yval-ycent).*50).^2)+(((xval-xcent).*200).^2)),'all','omitnan');
    
end
%group significance
orderedgroupsig=significantspread_groupsplit(E_MAP,:,:);
avgspreadwg=nan(9,500);

for group=1:9
    for itstimchn=1:size(orderedsigchns,2)
        tmp=reshape(squeeze(orderedgroupsig(1:64,itstimchn,group)),16,4);
        tmp=tmp(:,[1 3 4 2]);
        tmp(isnan(tmp))=0;
        xcent=sum(tmp.*[1 2 3 4],'all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
        ycent=sum(tmp.*(1:16)','all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
        
        xval=~isnan(tmp).*[1 2 3 4];
        xval(xval==0)=nan;
        yval=~isnan(tmp).*(1:16)';
        yval(yval==0)=nan;
        avgspreadwg(group,itstimchn)=mean(sqrt((((yval-ycent).*50).^2)+(((xval-xcent).*200).^2)),'all','omitnan');
        if ~isnan(xcent) && ~isnan(ycent)
            heatmap_centroid{ampit}(16-ycent+1:16-ycent+16,4-xcent+1:4-xcent+4,itstimchn+(group-1)*500)=tmp;
        end
    end
end

spreadgroup(ampit,1:9,1:length(avgspread))=avgspreadwg;
spread(ampit,1:length(avgspread))=avgspread;
figure(1000*ampit); hold on;cellfun(@(x) plot(-90:90,mean(x,1).*1000), groupdata{ampit})
color1 = linspace(0,1,groupit);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
indvsig=cellfun(@(x) size(x,1), groupdata{ampit});
text(-80,16,'# sig:')
text(-80,10,num2str(indvsig))
%text(-80,10,num2str(indvsig))V2
totalsig=sum(indvsig);
text(-80,20,['Total sig: ' num2str(totalsig) ' / ' num2str(totalchncount(ampit))])
xlabel('Time (ms)')
ylabel('Firing rate (Sp/s)')
xlim([-85 85])
%ylim([0 25])
ylim([0 400])
set(gca,'TickDir','out');
title([num2str(ampinterest) '\muA'])
end

%% Don't use
%%%%%%% binning based on 10uA trial only %doesn't work don't use
% trial avg single elect 
excitesupress=1;%1 for excite, 0 for supress - lower threshold for supress(see below)
cd([D_data(6).folder filesep D_data(6).name;])
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,[1,3,4,2]);
alldata=cell(length(savefilename{1}{2}.AMP),1);
trialsig=cell(9,length(savefilename{1}{2}.AMP),length(savefilename));
folderIDsig=cell(9,length(savefilename{1}{2}.AMP));
stimchncount=0;
spread=nan(length(savefilename{1}{2}.AMP),500);
spreadgroup=nan(length(savefilename{1}{2}.AMP),9,500);




ampinterest=[10 8];
groupdata=cell(8,1);
for i=1:8
    groupdata{i}=cell(9,1);
end
significantspread=nan(128,500);
totalchncount=0;
significantspread_groupsplit=nan(128,500,9);
for numerfolders=1:size(IDstructsave,1)
    trialend=length(fieldnames(IDstructsave{numerfolders}{1}));
    for trial=1:trialend
        ratespiking{numerfolders}(:,:,trial)=nanmean(IDstructsave{numerfolders}{3}{trial},3);
    end
   %%%%%%need to see if this line will pick the correct spiking trials.
   %%%%%%compare amp and select trial with that amp from savefile name.
   %%%%%%need to see if this owrrks with both the current steered data and
   %%%%%%non-current steering
    rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest(1)),1));

    totalchncount=64*size(rateAMPua,3)+totalchncount;
    meancatdata=cat(2,mean(rateAMPua(:,92:101,:),2),mean(rateAMPua(:,102:111,:),2),mean(rateAMPua(:,112:121,:),2),mean(rateAMPua(:,122:131,:),2),mean(rateAMPua(:,132:141,:),2),mean(rateAMPua(:,142:151,:),2),mean(rateAMPua(:,152:161,:),2),mean(rateAMPua(:,162:171,:),2),mean(rateAMPua(:,172:181,:),2));
    sig1=false(128,size(meancatdata,3),9);
    if excitesupress==0
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)-(std(rateAMPua(:,1:85,:),[],2).*1)); %%%%%NOTE LOWER THRESHOLD
        for groupit=1:9
            sig1(:,:,groupit)=squeeze(meancatdata(:,groupit,:))<thresh;
        end
        check_multiple=sum(sig1,3)>1;
        [~,pmax]=min(meancatdata,[],2);%pmin
        pmax=squeeze(pmax);%min
    else
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)+(std(rateAMPua(:,1:85,:),[],2).*3));
        for groupit=1:9
            sig1(:,:,groupit)=squeeze(meancatdata(:,groupit,:))>thresh;
        end
        check_multiple=sum(sig1,3)>1;
        [~,pmax]=max(meancatdata,[],2);
        pmax=squeeze(pmax);
    end
    

    rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest(2)),1));

    meancatdata=cat(2,mean(rateAMPua(:,92:101,:),2),mean(rateAMPua(:,102:111,:),2),mean(rateAMPua(:,112:121,:),2),mean(rateAMPua(:,122:131,:),2),mean(rateAMPua(:,132:141,:),2),mean(rateAMPua(:,142:151,:),2),mean(rateAMPua(:,152:161,:),2),mean(rateAMPua(:,162:171,:),2),mean(rateAMPua(:,172:181,:),2));
    sig2=false(128,size(meancatdata,3),9);
    if excitesupress==0
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)-(std(rateAMPua(:,1:85,:),[],2).*1)); %%%%%NOTE LOWER THRESHOLD
        for groupit=1:9
            sig2(:,:,groupit)=squeeze(meancatdata(:,groupit,:))<thresh;
        end
        check_multiple=sum(sig2,3)>1;
        [~,pmax]=min(meancatdata,[],2);%pmin
        pmax=squeeze(pmax);%min
    else
        thresh=squeeze(mean(rateAMPua(:,1:85,:),2)+(std(rateAMPua(:,1:85,:),[],2).*3));
        for groupit=1:9
            sig2(:,:,groupit)=squeeze(meancatdata(:,groupit,:))>thresh;
        end
        check_multiple=sum(sig2,3)>1;
        [~,pmax]=max(meancatdata,[],2);
        pmax=squeeze(pmax);
    end
  sig=(sig1 & sig2);
    
    for chn=1:64
        for stimchn=1:size(sig,2)
            if check_multiple(chn,stimchn)
                sig(chn,stimchn,:)=false;
                sig(chn,stimchn,pmax(chn,stimchn))=true;
            end
            group=find(squeeze(sig(chn,stimchn,:)==1));
            if ~isempty(group)
                for ampit=1:length(savefilename{1}{2}.AMP)
                    ampintt=savefilename{1}{2}.AMP(ampit);
                    rateselectedamp=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampintt),1));
                    groupdata{ampit}{group}=[groupdata{ampit}{group}; rateselectedamp(chn,:,stimchn)];
                    alldata{ampit}=[alldata{ampit}; rateselectedamp(chn,:,stimchn)];
                    significantspread(chn,stimchn+stimchncount)=meancatdata(chn,group,stimchn);
                    significantspread_groupsplit(chn,stimchn+stimchncount,group)=meancatdata(chn,group,stimchn);
                    trialsig{group,ampit,numerfolders}=[trialsig{group,ampit,numerfolders} stimchn];%determine how many different stim electrodes creating significance
                    folderIDsig{group,ampit}=[folderIDsig{group,ampit} str2double(savefilename{numerfolders}{3})];%determine how many different monkeys creating significance
                end
            end
            
        end
    end
    stimchncount=stimchncount+size(sig,2);
    
end
stimchncount=0;
%current significance
E_MAP=Depth;
orderedsigchns=significantspread(E_MAP,:);
avgspread=nan(size(orderedsigchns,2),1);

for itstimchn=1:size(orderedsigchns,2)
    tmp=reshape(orderedsigchns(1:64,itstimchn),16,4);
    tmp=tmp(:,[1 3 4 2]);
    tmp(isnan(tmp))=0;
    xcent=sum(tmp.*[1 2 3 4],'all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
    ycent=sum(tmp.*(1:16)','all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
 
    xval=~isnan(tmp).*[1 2 3 4];
    xval(xval==0)=nan;
    yval=~isnan(tmp).*(1:16)';
    yval(yval==0)=nan;
    avgspread(itstimchn)=mean(sqrt((((yval-ycent).*50).^2)+(((xval-xcent).*200).^2)),'all','omitnan');
    
end
%group significance
orderedgroupsig=significantspread_groupsplit(E_MAP,:,:);
avgspreadwg=nan(9,500);
for group=1:9
    for itstimchn=1:size(orderedsigchns,2)
        tmp=reshape(squeeze(orderedgroupsig(1:64,itstimchn,group)),16,4);
        tmp=tmp(:,[1 3 4 2]);
        tmp(isnan(tmp))=0;
        xcent=sum(tmp.*[1 2 3 4],'all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
        ycent=sum(tmp.*(1:16)','all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
        
        xval=~isnan(tmp).*[1 2 3 4];
        xval(xval==0)=nan;
        yval=~isnan(tmp).*(1:16)';
        yval(yval==0)=nan;
        avgspreadwg(group,itstimchn)=mean(sqrt((((yval-ycent).*50).^2)+(((xval-xcent).*200).^2)),'all','omitnan');
    end
end

spreadgroup(ampit,1:9,1:length(avgspread))=avgspreadwg;
spread(ampit,1:length(avgspread))=avgspread;
for ampit=1:5
figure(1000*ampit); hold on;cellfun(@(x) plot(-90:90,mean(x,1).*1000), groupdata{ampit})
color1 = linspace(0,1,groupit);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
indvsig=cellfun(@(x) size(x,1), groupdata{ampit});
text(-80,16,'# sig:')
text(-80,10,num2str(indvsig))
totalsig=sum(indvsig);
text(-80,20,['Total sig: ' num2str(totalsig) ' / ' num2str(totalchncount)])
xlabel('Time (ms)')
ylabel('Firing rate (Sp/s)')
xlim([-85 85])
ylim([0 25])
set(gca,'TickDir','out');
title([num2str(savefilename{1}{2}.AMP(ampit)) '\muA'])
end

%%
%%%%%%%%%%%%%%%%%% redundant %%%%%%%%%%%%%%%%%%%%% trial avg single elect
cd([D_data(6).folder filesep D_data(6).name;])
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,[1,3,4,2]);

ampinterest=10;
trialsig=cell(8,1);
folderIDsig=cell(8,1);
for iterate_time=10:10:80
ratespiking=cell(size(IDstructsave,1),1);
sigrate=[];
numstimelect=0;
folderIDsig{iterate_time/10}=[];%determine how many different monkeys creating significance
for numerfolders=1:size(IDstructsave,1) 
    trialend=length(fieldnames(IDstructsave{numerfolders}{1}));
     for trial=1:trialend
        ratespiking{numerfolders}(:,:,trial)=nanmean(IDstructsave{numerfolders}{3}{trial},3);
     end
    trialsig{iterate_time/10}{numerfolders}=[];%determine how many different stim electrodes creating significance
    rateAMPua=ratespiking{numerfolders}(:,:,find(savefilename{numerfolders}{2}.AMP==ampinterest):length(savefilename{numerfolders}{2}.AMP):end);
    numstimelect=size(rateAMPua,3)+numstimelect;
    for trial=1:size(rateAMPua,3)
        chnsignificant=zeros(16,4);
        for chn=1:64
            if mean(rateAMPua(chn,1:85,trial))+std(rateAMPua(chn,1:85,trial))*3<mean(rateAMPua(chn,92+(iterate_time-10):(92+iterate_time-1),trial))%ttest(rateAMPua(chn,randperm(85,10),trial),rateAMPua(chn,92+(iterate_time-10):92+(iterate_time-1),trial),'Tail','left')==1%%&&any(movmean(rateAMPua(chn,:,trial),5).*1000>20)%),"Tail","left")==1%
                sigrate=cat(1,sigrate,rateAMPua(chn,:,trial));
                chnsignificant(ordershapearray==chn)=1;
%                 figure (1)
%                 hold on
%                 plot(-90:90,movmean(rateAMPua(chn,:,trial),5).*1000)
%                 x=0;
            trialsig{iterate_time/10}{numerfolders}=[trialsig{iterate_time/10}{numerfolders} trial];%determine how many different stim electrodes creating significance
            folderIDsig{iterate_time/10}=[folderIDsig{iterate_time/10} str2double(savefilename{numerfolders}{3})];%determine how many different monkeys creating significance
            end
        end
%         if sum(chnsignificant,'all')>16
%         figure
%         heatmap(chnsignificant);
%         end
    end
end
text(-80,13-(iterate_time/10),num2str(size(sigrate,1)))
size(sigrate,1)
figure(5000); plot(-90:90,mean(sigrate).*1000)
hold on
end
text(-80,14,'# sig:')
xlabel('Time (ms)')
ylabel('Firing rate (Sp/s)')
xlim([-85 85])
color1 = linspace(0,1,length(10:10:80));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
set(gca,'TickDir','out');
%%
%%%%%%%%%%%%%%%%%% redundant %%%%%%%%%%%%%%%%%%%%%%%% # Sig monkeys and stim elect per time
sigmonkeys=cellfun(@(x) unique(x), folderIDsig, 'UniformOutput', false);
sigmonkeys=cellfun(@(x) sum(diff(x)>3)+1,sigmonkeys);

for timeloop=1:8
sigstimelect{timeloop,1}=sum(cellfun(@(x) length(unique(x)),trialsig{timeloop},'UniformOutput',true));
end
%% sort into responding and not responding channels 
ampinterest=6;
ratespiking_nill=cell(size(storefileamp,1),1);
ratespiking_excite=cell(size(storefileamp,1),1);
ratespiking_supress=cell(size(storefileamp,1),1);
for numfolders=1:size(storefileamp,1)
    ampit=reshape(storefileamp(numfolders,:),[],size(storefilechn,2));%%%%%%check this!!! need to see if moving has affected
    if any(ampit~=-1 & ampit~=ampinterest,'all')% pull out 6uA results
        continue
    end
    for chn=65:128
        tmp=squeeze(ratespiking{numfolders}(chn,:,255));
        %         epoch=10;%ms
        %         TimeWindow=0;
        %         for time=0:floor(90/epoch)-2
        %             if  ttest(tmp(75:85),tmp(92+(time*epoch):92+((time+1)*epoch)))==1
        %                 TimeWindow=92+(time*epoch):92+((time+1)*epoch);
        %                 break
        %             end
        %         end
        
        if ttest(tmp(1:85),tmp(92:92+84),"Tail","left")==1
            ratespiking_excite{numfolders}=cat(1,ratespiking_excite{numfolders}, ratespiking{numfolders}(chn,:,:));
        elseif ttest(tmp(1:85),tmp(92:92+84),"Tail","right")==1
            ratespiking_supress{numfolders}=cat(1,ratespiking_supress{numfolders}, ratespiking{numfolders}(chn,:,:));
        else
            ratespiking_nill{numfolders}=cat(1,ratespiking_nill{numfolders}, ratespiking{numfolders}(chn,:,:));
        end
    end
end
%% sort into responding not responding using standard deviations
ampinterest=6;
ratespiking_nill=cell(size(storefileamp,1),1);
ratespiking_excite=cell(size(storefileamp,1),1);
ratespiking_supress=cell(size(storefileamp,1),1);
for numfolders=1:size(storefileamp,1)
    ampit=reshape(storefileamp(numfolders,:),[],size(storefilechn,2));%%%%%%check this!!! need to see if moving has affected
    if any(ampit~=-1 & ampit~=ampinterest,'all')% pull out 6uA results
        continue
    end
    for chn=1:64
        tmp=squeeze(ratespiking{numfolders}(chn,:,255));
        %         epoch=10;%ms
        %         TimeWindow=0;
        %         for time=0:floor(90/epoch)-2
        %             if  ttest(tmp(75:85),tmp(92+(time*epoch):92+((time+1)*epoch)))==1
        %                 TimeWindow=92+(time*epoch):92+((time+1)*epoch);
        %                 break
        %             end
        %         end
        
        if any((mean(tmp(1:85))+std(tmp(1:85))*3)<tmp(92+0:92+8))%92:92+84
            ratespiking_excite{numfolders}=cat(1,ratespiking_excite{numfolders}, ratespiking{numfolders}(chn,:,:));
        elseif any((mean(tmp(1:85))-std(tmp(1:85)))>tmp(92:92+84))
            ratespiking_supress{numfolders}=cat(1,ratespiking_supress{numfolders}, ratespiking{numfolders}(chn,:,:));
        else
            ratespiking_nill{numfolders}=cat(1,ratespiking_nill{numfolders}, ratespiking{numfolders}(chn,:,:));
        end
    end
end


%% sort into num elects
ratespk=ratespiking_excite;
spikingperactivechn=cell(size(storefilechn,2),1);
rateperactivechn=cell(size(storefilechn,2),1);
for iteratefile=1:size(storefileamp,1)
    if ~isempty(ratespk{iteratefile})
        for iteratechnsamp=1:size(ampit,1)
            numchnsactive=sum(ampit(iteratechnsamp,:)==ampinterest);
            spikingperactivechn{numchnsactive}=[spikingperactivechn{numchnsactive} avgspiking{iteratefile}(:,iteratechnsamp)];
            rateperactivechn{numchnsactive}=cat(1,rateperactivechn{numchnsactive},ratespk{iteratefile}(:,:,iteratechnsamp));
        end
    end
end

%plotting
rateplot_collapsedchn=cellfun(@(x) nanmean(x,1),rateperactivechn,'UniformOutput',false);
figure
hold on
title('V2 rate @ 6\muA')
cellfun(@(x) plot(-90:90,x.*1000),rateplot_collapsedchn','UniformOutput',false)
color1 = linspace(0,1,length(rateplot_collapsedchn));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
lgd = legend('1','2','3','4','5','6','7','8');
lgd.Title.String = '# stim elect';
ylabel('Firing rate (Sp/s)')
xlabel('Time (ms)')
set(gca,'TickDir','out');
%%
%% sort into which channel is responding 
% sort into responding and not responding channels
chn_interest=1:64;
excitesupnill=1;%excite=1,suppress=-1,nill=0
ampinterest=6;
ratespk=ratespiking;
spikingperactivechn=cell(size(storefilechn,2),1);
rateperactivechn=cell(size(storefilechn,2),1);
loop=0;
for iteratefile=3%:size(storefileamp,1)
    ampit=reshape(storefileamp(iteratefile,:),[],size(storefilechn,2));%%%%%%check this!!! need to see if moving has affected
    if any(ampit~=-1 & ampit~=ampinterest,'all')% pull out 6uA results
        continue
    end
    if ~isempty(ratespk{iteratefile})
        for iteratechnsamp=1:8%size(ampit,1)
            numchnsactive=sum(ampit(iteratechnsamp,:)==ampinterest);
            spikingperactivechn{numchnsactive}=[spikingperactivechn{numchnsactive} avgspiking{iteratefile}(:,iteratechnsamp)];
            for iteratechn=chn_interest

                tmp=ratespk{iteratefile}(iteratechn,:,iteratechnsamp);
                peakpos=find(tmp(92:end)==max(tmp(92:end)),1,'first')-10:find(tmp(92:end)==max(tmp(92:end)),1,'first')+10;
                if any(peakpos<1)
                    peakpos=1:21;
                elseif any(peakpos+91>181)
                    peakpos=65:85;
                end
%                         epoch=10;%ms
%                         TimeWindow=0;
%                         for time=0:floor(90/epoch)-2
%                             if  ttest(tmp(75:85),tmp(92+(time*epoch):92+((time+1)*epoch)))==1
%                                 TimeWindow=92+(time*epoch):92+((time+1)*epoch);
%                                 break
%                             end
%                         end

                if excitesupnill==1 &&  (ttest(tmp(65:85),tmp(peakpos+91),"Tail","left")==1)% || ttest(tmp(65:85),tmp(23:43),"Tail","left")==1 || ttest(tmp(65:85),tmp(44:64),"Tail","left")==1 || ttest(tmp(65:85),tmp(65:85),"Tail","left")==1)
                    loop=loop+1;
                    rateperactivechn{numchnsactive}=cat(1,rateperactivechn{numchnsactive},tmp);
                elseif excitesupnill==-1 &&  (ttest(tmp(65:85),tmp(2:22),"Tail","right")==1 || ttest(tmp(65:85),tmp(23:43),"Tail","right")==1 || ttest(tmp(65:85),tmp(44:64),"Tail","right")==1 || ttest(tmp(65:85),tmp(65:85),"Tail","right")==1)
                    rateperactivechn{numchnsactive}=cat(1,rateperactivechn{numchnsactive},tmp);
                elseif excitesupnill==0 && ttest(tmp(1:85),tmp(92:92+84))==0
                    rateperactivechn{numchnsactive}=cat(1,rateperactivechn{numchnsactive},tmp);
                end
            end
        end
    end
end
%plotting
rateplot_collapsedchn=cellfun(@(x) nanmean(x,1),rateperactivechn,'UniformOutput',false);
figure
hold on
title('V2 rate @ 6\muA')
cellfun(@(x) plot(-90:90,x.*1000),rateplot_collapsedchn','UniformOutput',false)
color1 = linspace(0,1,length(rateplot_collapsedchn));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
lgd = legend('1','2','3','4','5','6','7','8');
lgd.Title.String = '# stim elect';
ylabel('Firing rate (Sp/s)')
xlabel('Time (ms)')
set(gca,'TickDir','out');
%% sort based on each trial not just all 8 elect on
% sort into responding and not responding channels
chn_interest=65:128;
excitesupnill=1;%excite=1,suppress=-1,nill=0
ampinterest=6;
ratespk=ratespiking;
spikingperactivechn=cell(size(storefilechn,2),1);
rateperactivechn=cell(size(storefilechn,2),1);
for iteratefile=1:size(storefileamp,1)
    ampit=reshape(storefileamp(iteratefile,:),[],size(storefilechn,2));%%%%%%check this!!! need to see if moving has affected
    if any(ampit~=-1 & ampit~=ampinterest,'all')% pull out 6uA results
        continue
    end
    if ~isempty(ratespk{iteratefile})
        for iteratechnsamp=1:size(ampit,1)
            numchnsactive=sum(ampit(iteratechnsamp,:)==ampinterest);
            spikingperactivechn{numchnsactive}=[spikingperactivechn{numchnsactive} avgspiking{iteratefile}(:,iteratechnsamp)];
            for iteratechn=chn_interest

                tmp=ratespk{iteratefile}(iteratechn,:,iteratechnsamp);
%                         epoch=10;%ms
%                         TimeWindow=0;
%                         for time=0:floor(90/epoch)-2
%                             if  ttest(tmp(75:85),tmp(92+(time*epoch):92+((time+1)*epoch)))==1
%                                 TimeWindow=92+(time*epoch):92+((time+1)*epoch);
%                                 break
%                             end
%                         end
%%%%%% check the conditions
%%%%%% below!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if excitesupnill==1 &&  (ttest(tmp(1:85),tmp(92:92+84),"Tail","left")==1)% || ttest(tmp(65:85),tmp(23:43),"Tail","left")==1 || ttest(tmp(65:85),tmp(44:64),"Tail","left")==1 || ttest(tmp(65:85),tmp(65:85),"Tail","left")==1)
                    rateperactivechn{numchnsactive}=cat(1,rateperactivechn{numchnsactive},tmp);
                elseif excitesupnill==-1 &&  (ttest(tmp(65:85),tmp(2:22),"Tail","right")==1 || ttest(tmp(65:85),tmp(23:43),"Tail","right")==1 || ttest(tmp(65:85),tmp(44:64),"Tail","right")==1 || ttest(tmp(65:85),tmp(65:85),"Tail","right")==1)
                    rateperactivechn{numchnsactive}=cat(1,rateperactivechn{numchnsactive},tmp);
                elseif excitesupnill==0 && ttest(tmp(1:85),tmp(92:92+84))==0
                    rateperactivechn{numchnsactive}=cat(1,rateperactivechn{numchnsactive},tmp);
                end
            end
        end
    end
end
%plotting
rateplot_collapsedchn=cellfun(@(x) nanmean(x,1),rateperactivechn,'UniformOutput',false);
figure
hold on
title('V2 rate @ 6\muA')
cellfun(@(x) plot(-90:90,x.*1000),rateplot_collapsedchn','UniformOutput',false)
color1 = linspace(0,1,length(rateplot_collapsedchn));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
lgd = legend('1','2','3','4','5','6','7','8');
lgd.Title.String = '# stim elect';
ylabel('Firing rate (Sp/s)')
xlabel('Time (ms)')
set(gca,'TickDir','out');
%%

%% plotting
spikingperelect=cellfun(@(x) nanmean(x(65:128,:),'all'),spikingperactivechn);
figure
plot(spikingperelect.*1000)
title('V1 @ 6\muA')
ylabel('Firing rate (Sp/s)')
xlabel('# stim elect')

spikingperelect=cellfun(@(x) nanmean(x(1:64,:),'all'),spikingperactivechn);
figure
plot(spikingperelect.*1000)
title('V2 @ 6\muA')
ylabel('Firing rate (Sp/s)')
xlabel('# stim elect')
%rateplot=cellfun(@(x) nanmean(x(65:128,:,:),1),rateperactivechn,'UniformOutput',false);
rateplot_collapsedchn=cellfun(@(x) nanmean(x,1),rateperactivechn,'UniformOutput',false);
figure
hold on
title('V1 rate @ 6\muA')
cellfun(@(x) plot(-90:90,x.*1000),rateplot_collapsedchn','UniformOutput',false)
color1 = linspace(0,1,length(rateplot_collapsedchn));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
lgd = legend('1','2','3','4','5','6','7','8');
lgd.Title.String = '# stim elect';
ylabel('Firing rate (Sp/s)')
xlabel('Time (ms)')

rateplot=cellfun(@(x) nanmean(x(1:64,:,:),1),rateperactivechn,'UniformOutput',false);
rateplot_collapsedchn=cellfun(@(x) mean(x,1),rateplot,'UniformOutput',false);

figure
hold on
title('V2 rate @ 6\muA')
cellfun(@(x) plot(-90:90,x.*1000),rateplot_collapsedchn','UniformOutput',false)
color1 = linspace(0,1,length(rateplot_collapsedchn));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
lgd = legend('1','2','3','4','5','6','7','8');
lgd.Title.String = '# stim elect';
ylabel('Firing rate (Sp/s)')
xlabel('Time (ms)')
%% Which folders are weird?
rateplot=cellfun(@(x) nanmean(x(1:64,:,:),3),ratespiking,'UniformOutput',false);
figure
hold on
title('V1 rate') 
cellfun(@(x) plot(nanmean(x,1).*1000),rateplot','UniformOutput',false)
lgd = legend('1','2','3','4','5','6','7','8');
lgd.Title.String = 'folder';
color1 = linspace(0,1,length(rateplot_collapsedchn));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
lgd = legend('1','2','3','4','5','6','7','8');
lgd.Title.String = 'folder';
%folder 4, 3 and 2 are bad