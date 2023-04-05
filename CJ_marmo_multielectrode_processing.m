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

        [IDstruct, baslinespikestruct,ratestruct] = sortTrials_SM(timestart,timeend,trig,0,1,1,endtrial,chn_range);
        stimfilename=dir('*exp_datafile_*');
        IDstructsave{k}={IDstruct baslinespikestruct ratestruct};
        stimVar=load(stimfilename.name,'AMP','CHN');
        savefilename{k}=[{str2double(stimfilename.name(end-6:end-4))} {stimVar}];

end
%folders=D_data(checkfolder,:);

IDstructsave=IDstructsave(checkfolder,:);
savefilename=savefilename(checkfolder);
%% Concat trials with same stimulus file
%filenameunique=unique(savefilename);
checkoff=false(length(savefilename),1);
IDstructcompiled=cell(length(IDstructsave),1);
BaselineIDstructcompiled=cell(length(IDstructsave),1);
ratesructcompiled=cell(length(IDstructsave),1);
storefilechn=zeros(length(IDstructsave),8);
storefileamp=zeros(length(IDstructsave),2040);
uniqueamp=0;
for loop=1:length(IDstructsave)% loop through the stimulation pairs. Avoid using the first ones
     trialend=length(fieldnames(IDstructsave{loop}{1}));
    filessame=all(storefilechn==savefilename{loop}{2}.CHN,2) & all(storefileamp==savefilename{loop}{2}.AMP(:)',2);
    if ~any(filessame)
        uniqueamp=uniqueamp+1;
        storefileamp(uniqueamp,:)=savefilename{loop}{2}.AMP(:)';
        storefilechn(uniqueamp,:)=savefilename{loop}{2}.CHN;
        for trial=1:trialend
            tnum=['T', num2str(trial)];
            if ~isempty(IDstructsave{loop}{1}.(tnum))
             numnantoadd=128-size(IDstructsave{loop}{1}.(tnum),1);
            ids_full=[IDstructsave{loop}{1}.(tnum); nan(numnantoadd,size(IDstructsave{loop}{1}.(tnum),2))];
            bss_full=[IDstructsave{loop}{2}.(tnum); nan(numnantoadd,size(IDstructsave{loop}{2}.(tnum),2))];
            rate_full=[IDstructsave{loop}{3}{trial};nan(numnantoadd,size(IDstructsave{loop}{3}{trial},2),size(IDstructsave{loop}{3}{trial},3))];
            IDstructcompiled{uniqueamp}.(tnum)=ids_full;
            BaselineIDstructcompiled{uniqueamp}.(tnum)=bss_full;
            ratesructcompiled{uniqueamp}.(tnum)=rate_full;
            end
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
            ratesructcompiled{index}.(tnum)=cat(3,ratesructcompiled{index}.(tnum), rate_full);
            end
        end
    end
end


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

%% sort into num elects
ampinterest=6;
spikingperactivechn=cell(size(storefilechn,2),1);
rateperactivechn=cell(size(storefilechn,2),1);
for iteratefile=1:size(storefileamp,1)
    ampit=reshape(storefileamp(iteratefile,:),[],size(storefilechn,2));
    if any(ampit~=-1 & ampit~=ampinterest,'all')% pull out 6uA results
        continue
    end
    for iteratechnsamp=1:size(ampit,1)
        numchnsactive=sum(ampit(iteratechnsamp,:)==ampinterest);
        spikingperactivechn{numchnsactive}=[spikingperactivechn{numchnsactive} avgspiking{iteratefile}(:,iteratechnsamp)];
        rateperactivechn{numchnsactive}=cat(3,rateperactivechn{numchnsactive},ratespiking{iteratefile}(:,:,iteratechnsamp));
    end
end
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
rateplot=cellfun(@(x) nanmean(x(65:128,:,:),3),rateperactivechn,'UniformOutput',false);
rateplot_collapsedchn=cellfun(@(x) nanmean(x,1),rateplot,'UniformOutput',false);
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

rateplot=cellfun(@(x) nanmean(x(1:64,:,:),3),rateperactivechn,'UniformOutput',false);
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
rateplot=cellfun(@(x) nanmean(x(65:128,:,:),3),ratespiking,'UniformOutput',false);
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