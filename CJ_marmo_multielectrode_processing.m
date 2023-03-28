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
try
       [IDstruct, baslinespikestruct] = sortTrials_SM(timestart,timeend,trig,0,1,1,endtrial);
       stimfilename=dir('*exp_datafile_*');
       IDstructsave{k}={IDstruct baslinespikestruct};
       stimVar=load(stimfilename.name,'AMP','CHN');
       savefilename{k}=[{str2double(stimfilename.name(end-6:end-4))} {stimVar}];
catch
    continue
end
end
%folders=D_data(checkfolder,:);

IDstructsave=IDstructsave(checkfolder,:);
savefilename=savefilename(checkfolder);
%% Concat trials with same stimulus file
%filenameunique=unique(savefilename);
checkoff=false(length(savefilename),1);
IDstructcompiled=cell(length(IDstructsave),1);
BaselineIDstructcompiled=cell(length(IDstructsave),1);
loop=0;
     
     storefilechn=zeros(length(IDstructsave),8);
    storefileamp=zeros(length(IDstructsave),2040);
    uniqueamp=0;
for loop=1:length(IDstructsave)% loop through the stimulation pairs. Avoid using the first ones
%    currD = folders(loop).name; % Get the current subdirectory name
%     cd([folders(loop).folder filesep currD])
    
    
    
    %k=find(checkoff==false,1);
    filessame=all(storefilechn==savefilename{loop}{2}.CHN,2) & all(storefileamp==savefilename{loop}{2}.AMP(:)',2);
    if ~any(filessame)
        uniqueamp=uniqueamp+1;
        storefileamp(uniqueamp,:)=savefilename{loop}{2}.AMP(:)';
        storefilechn(uniqueamp,:)=savefilename{loop}{2}.CHN;
        for trial=1:trialend
            tnum=['T', num2str(trial)];
            
            IDstructcompiled{uniqueamp}.(tnum)=IDstructsave{loop}{1}.(tnum);
            BaselineIDstructcompiled{uniqueamp}.(tnum)=IDstructsave{loop}{2}.(tnum);
            
        end
    else
        index=find(filessame);
        for trial=1:trialend
            tnum=['T', num2str(trial)];
        IDstructcompiled{index}.(tnum)=[IDstructcompiled{index}.(tnum) IDstructsave{loop}{1}.(tnum)];
        BaselineIDstructcompiled{index}.(tnum)=[BaselineIDstructcompiled{index}.(tnum) IDstructsave{loop}{2}.(tnum)];
        end
    end
%     filessame=any(savefilename==filenameunique(loop),2);
%     checkoff=checkoff | filessame;
%     %[avgnospT,stderrspktrial,~] = AverageTrialResponse_SM(IDstruct, baslinespikestruct);
%     index=find(filessame==true);
%     trialend=length(fieldnames(IDstructsave{index(1)}{1}));
%     for loopfiles=1:length(index)
%         for trial=1:trialend
%             tnum=['T', num2str(trial)];
%             if loopfiles==1
%                 IDstructcompiled{loop}.(tnum)=IDstructsave{index(loopfiles)}{1}.(tnum);
%                 BaselineIDstructcompiled{loop}.(tnum)=IDstructsave{index(loopfiles)}{2}.(tnum);
%             else
%                 IDstructcompiled{loop}.(tnum)=[IDstructcompiled{loop}.(tnum) IDstructsave{index(loopfiles)}{1}.(tnum)];
%                 BaselineIDstructcompiled{loop}.(tnum)=[BaselineIDstructcompiled{loop}.(tnum) IDstructsave{index(loopfiles)}{2}.(tnum)];
%             end
%         end
%     end
end


%% trial average
avgspiking=cell(loop,1);
for numerfolders=1:loop
    [avgnospT,stderrspktrial,~]=AverageTrialResponse_SM(IDstructcompiled{numerfolders},BaselineIDstructcompiled{numerfolders});
    avgspiking{numerfolders}=avgnospT;
end

%% bin into #stim elect

