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
%% trial avg single elect


ampinterest=10;
cd([D_data(6).folder filesep D_data(6).name;])
ratespiking=cell(size(storefileamp,1),1);
sigrate=[];
for numerfolders=1:size(storefileamp,1)
     for trial=1:trialend
        ratespiking{numerfolders}(:,:,trial)=nanmean(IDstructsave{numerfolders}{3}{trial},3);
     end
    rateAMPua=ratespiking{numerfolders}(:,:,find(savefilename{numerfolders}{2}.AMP==ampinterest):length(savefilename{numerfolders}{2}.AMP):end);
    for trial=1:size(rateAMPua,3)
        for chn=1:64
            if any(std(rateAMPua(chn,1:85,trial))*15<(rateAMPua(chn,92+10:92+30,trial)))&&any(movmean(rateAMPua(chn,:,trial),5).*1000>20)%),"Tail","left")==1%
                sigrate=cat(1,sigrate,rateAMPua(chn,:,trial));
                figure (1)
                hold on
                plot(-90:90,movmean(rateAMPua(chn,:,trial),5).*1000)
                x=0;
            end
        end
    end
end


%% sort into responding and not responding channels - ttest
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