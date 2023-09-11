function V2respcomb=multielectplot(ratestruct,savefilename,chnrange,excitesupress, chnexclude)
ampinterest=6;

 s = RandStream('mt19937ar','Seed',296);
% pre-processing
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=cell(numfolderstotal,1);

for numerfolders=1:numfolderstotal
    trialend=size(savefilename{numerfolders}{2}.AMP,1);
    for trial=1:trialend
        trialcheck=['T', num2str(trial)];
        ratespiking{numerfolders}(:,:,trial)=nanmean(ratestruct{numerfolders}.(trialcheck),3);
        %make sure only responding elect included
        subselectbaseline=1:89;
        permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
        subselectbaseline=permuted_numbers(1:50);
        thresh=squeeze(mean(ratespiking{numerfolders}(:,subselectbaseline,trial),2)+(std(ratespiking{numerfolders}(:,subselectbaseline,trial),[],2).*3));
        notsigsingle=squeeze(mean(ratespiking{numerfolders}(:,92:181,trial),2))<=thresh;
        ratespiking{numerfolders}(notsigsingle,:,trial)=nan;
    end
end

%%
% 1. # elect vs FR.
numelecttotal=length(savefilename{1}{2}.CHN);
separateeachstimelect=cell(numelecttotal,numelecttotal);% 8 plots with 8 different #elect
respelect=cell(numelecttotal,1);
V2resp=cell(numelecttotal,1);
V1resp=cell(numelecttotal,1);
for numstimelect=1:numelecttotal
    for numerfolders=1:numfolderstotal
        trialend=size(savefilename{numerfolders}{2}.AMP,1);
        for trial=1:trialend
            whichelectstim=savefilename{numerfolders}{2}.AMP(trial,:)==ampinterest;
            if sum(whichelectstim)~=numstimelect%check if the trial has the correct num stim elect and correct amplitude stim
                continue
            end
            index=find(whichelectstim);
          for i=1:numstimelect
            separateeachstimelect{index(i),numstimelect}=cat(1,separateeachstimelect{index(i),numstimelect},ratespiking{numerfolders}(1:64,:,trial));
          end
            diff(squeeze(mean(ratespiking{numerfolders}(:,92:181,:),2)));
            V2resp{numstimelect}=cat(1,V2resp{numstimelect},ratespiking{numerfolders}(1:64,:,trial));
            V1resp{numstimelect}=cat(1,V1resp{numstimelect},ratespiking{numerfolders}(65:128,:,trial));
            respelect{numstimelect}=cat(3,respelect{numstimelect},ratespiking{numerfolders}(:,:,trial));
        end
    end
end
%% sigmoid

timetolook=92:181;
timetolookV1=92:112;
s = RandStream('mt19937ar','Seed',296);
subselectbaseline=1:89;
permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
subselectbaseline=permuted_numbers(1:50);


%V2dat=cellfun(@(x) mean(max(x(:,timetolook)-mean(x(:,subselectbaseline),2,'omitnan'),[],2),'omitnan'),V2resp).*1000;
V2dat=cellfun(@(x) mean(x(:,timetolook)-mean(x(:,subselectbaseline),2,'omitnan'),'all','omitnan'),V2resp).*1000;
sdV2resp=cellfun(@(x) SEM(x(:,timetolook)-mean(x(:,subselectbaseline),2,'omitnan'),0),V2resp).*1000;
V1dat=cellfun(@(x) mean(x(:,timetolookV1)-mean(x(:,subselectbaseline),2,'omitnan'),'all','omitnan'),V1resp).*1000;
sdV1resp=cellfun(@(x) SEM(x(:,timetolookV1)-mean(x(:,subselectbaseline),2,'omitnan'),0),V1resp).*1000;
figure
hold on
errorbar(1:8,V2dat,sdV2resp,'b')
set(gca,'TickDir','out');
title('V2')
xlabel('Number of stim electrodes in V1')
ylabel('Firing rate (sp/s)')
beautifyPlot;

figure
hold on
errorbar(1:8,V1dat,sdV1resp,'r')
set(gca,'TickDir','out');
title('V1')
xlabel('Number of stim electrodes in V1')
ylabel('Sp/s')
beautifyPlot;
%significance
 V2data=cellfun(@(x) mean(x(:,timetolook)-mean(x(:,subselectbaseline),2,'omitnan'),2,'omitnan').*1000,V2resp,'UniformOutput',false);
dat=nan(8960,8);
for num=1:8
dat(1:length(V2data{num}),num)=V2data{num};
end
%anova1(dat)
%num active elect

numelectpertrial=cellfun(@(x) mean(squeeze(sum(mean(x(1:64,timetolook,:),2,'omitnan')>0,1,'omitnan'))),respelect,'UniformOutput',true);
numelectpertrialSD=cellfun(@(x) SEM(squeeze(sum(mean(x(1:64,timetolook,:),2,'omitnan')>0,1,'omitnan')),0),respelect,'UniformOutput',true);
figure
errorbar(1:8,numelectpertrial,numelectpertrialSD)
title('V2')
ylabel('Number active electrodes in V2')
xlabel('Number of stim electrodes in V1')
set(gca,'TickDir','out');
beautifyPlot

%sigmoide combines groups 1-2, 3-4, 5-6 7-8
V2respcomb=cell(4,1);
V2respelect=cell(4,1);
V1respcomb=cell(4,1);

for i=1:length(V2resp)/2
V2respcomb{i}=cat(1,V2resp{i*2-1},V2resp{i*2});
V2respelect{i}=cat(3,respelect{i*2-1},respelect{i*2});
V1respcomb{i}=cat(1,V1resp{i*2-1},V1resp{i*2});

end
V2datcomb=cellfun(@(x) mean(x(:,timetolook)-mean(x(:,subselectbaseline),2,'omitnan'),'all','omitnan'),V2respcomb).*1000;
V1datcomb=cellfun(@(x) mean(x(:,timetolookV1)-mean(x(:,subselectbaseline),2,'omitnan'),'all','omitnan'),V1respcomb).*1000;

testsig=cellfun(@(x) mean(x(:,timetolook)-mean(x(:,subselectbaseline),2,'omitnan'),2,'omitnan'),V2respcomb,'UniformOutput',false);
testsigv2=nan(17000,4);
for i=1:4
testsigv2(1:length(testsig{i}(testsig{i}>0)),i)=testsig{i}(testsig{i}>0);
end
sdV2respcomb=cellfun(@(x) SEM(x(:,timetolook)-mean(x(:,subselectbaseline),2,'omitnan'),0),V2respcomb).*1000;
sdV1respcomb=cellfun(@(x) SEM(x(:,timetolookV1)-mean(x(:,subselectbaseline),2,'omitnan'),0),V1respcomb).*1000;
figure
hold on
errorbar(1:2:8,V2datcomb,sdV2respcomb,'b')
ylim([0.5 1.1])
set(gca,'TickDir','out');
title('V2')
xlabel('Number of stim electrodes in V1')
ylabel('Firing rate (sp/s)')
beautifyPlot;
 cellfun(@(x) sum(~isnan(x(:,1))),V2respcomb)
 
 figure
hold on
errorbar(1:2:8,V1datcomb,sdV1respcomb,'b')
set(gca,'TickDir','out');
title('V1')
xlabel('Number of stim electrodes in V1')
ylabel('Firing rate (sp/s)')
beautifyPlot;

  cellfun(@(x) sum(~isnan(x(:,1))),V1respcomb)
%numactive elect
numelectpertrial=cellfun(@(x) mean(squeeze(sum(mean(x(1:64,timetolook,:),2,'omitnan')>0,1,'omitnan'))),V2respelect,'UniformOutput',true);
numelectpertrialSD=cellfun(@(x) SEM(squeeze(sum(mean(x(1:64,timetolook,:),2,'omitnan')>0,1,'omitnan')),0),V2respelect,'UniformOutput',true);
figure
errorbar(1:2:8,numelectpertrial,numelectpertrialSD)
title('V2')
ylabel('Number active electrodes in V2')
xlabel('Number of stim electrodes in V1')
set(gca,'TickDir','out');
beautifyPlot


%% compare single and dual electrode centroid position
areaV=2;
if areaV==2
    chns=1:64;
    cols=1:4;
else
    chns=65:128;
    cols=5:8;
end
timetolook=92:181;
stimarray=[1:16;33:48;49:64;17:32;]'+64;
s = RandStream('mt19937ar','Seed',296);
subselectbaseline=1:89;
permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
subselectbaseline=permuted_numbers(1:50);

E_MAP=Depth(1);
emap_shaped=reshape(E_MAP,16,8);
emap_shaped=emap_shaped(:,[1 3 4 2 5 7 8 6]);

V2resp_centroid=cell(6,1);
numstimelect=2;

for numerfolders=1:numfolderstotal
    trialend=size(savefilename{numerfolders}{2}.AMP,1);
    for trial=1:trialend
        whichelectstim=savefilename{numerfolders}{2}.AMP(trial,:)==ampinterest;
        if sum(whichelectstim)~=numstimelect%check if the trial has the correct num stim elect and correct amplitude stim
            continue
        end
        
        E1=find(whichelectstim,1,'first');
        E2=find(whichelectstim,1,'last');
        trialE1=find(savefilename{numerfolders}{2}.AMP(:,E1)==ampinterest,1,'first'); 
        trialE2=find(savefilename{numerfolders}{2}.AMP(:,E2)==ampinterest,1,'first');
        %V2dual=max((ratespiking{numerfolders}(1:64,timetolook,trial)-mean(ratespiking{numerfolders}(1:64,subselectbaseline,trial),2,'omitnan')),[],2,'omitnan');
        V2dual=mean(ratespiking{numerfolders}(chns,timetolook,trial)-mean(ratespiking{numerfolders}(chns,subselectbaseline,trial),2,'omitnan'),2,'omitnan');
        V2E1=mean(ratespiking{numerfolders}(chns,timetolook,trialE1)-mean(ratespiking{numerfolders}(chns,subselectbaseline,trial),2,'omitnan'),2,'omitnan');
        V2E2=mean(ratespiking{numerfolders}(chns,timetolook,trialE2)-mean(ratespiking{numerfolders}(chns,subselectbaseline,trial),2,'omitnan'),2,'omitnan');


            [centV2rdual,centV2cdual]=centroidpos(V2dual(emap_shaped(:,1:4)));
            [centV2rE1,centV2cE1]=centroidpos(V2E1(emap_shaped(:,1:4)));
            [centV2rE2,centV2cE2]=centroidpos(V2E2(emap_shaped(:,1:4)));
%              [centV2rdual,centV2cdual]=find(V2dual(emap_shaped(:,1:4))==max(V2dual),1,'first');
%             [centV2rE1,centV2cE1]=find(V2E1(emap_shaped(:,1:4))==max(V2E1),1,'first');
%             [centV2rE2,centV2cE2]=find(V2E2(emap_shaped(:,1:4))==max(V2E2),1,'first');
            rdist=(abs(centV2rdual-centV2rE1)+abs(centV2rdual-centV2rE2))/2;
            cdist=(abs(centV2cdual-centV2cE1)+abs(centV2cdual-centV2cE2))/2;
%             if isempty(rdist) || isempty(cdist) 
%                 cdist=nan;
%                 rdist=nan;
%             end
            V2resp_centroid{6}=cat(1,V2resp_centroid{6},[rdist,cdist]);
            E1num=savefilename{numerfolders}{1, 4}(trialE1,3);  
             E2num=savefilename{numerfolders}{1, 4}(trialE2,3);
             [r1,c1]=find(stimarray==E1num);
             [r2,c2]=find(stimarray==E2num);
             rdistV1=abs(r1-r2);
             cdistV1=abs(c1-c2);
            V2resp_centroid{2}=cat(1,V2resp_centroid{2},[rdistV1,cdistV1]);
            %for num elect
            V2resp_centroid{1}=cat(1,V2resp_centroid{1},sum(V2dual>0));%%
            %for num elect single
            V2resp_centroid{4}=cat(1,V2resp_centroid{4},(sum(V2E1>0)+sum(V2E2>0)));%%
            %for spiking rate
            V2resp_centroid{3}=cat(1,V2resp_centroid{3},mean(V2dual(V2dual>0),'omitnan'));
            %for spiking rate single
            V2resp_centroid{5}=cat(1,V2resp_centroid{5},mean(V2E1(V2E1>0),'omitnan')+mean(V2E2(V2E2>0),'omitnan'));
    end
end

ElectdistV1=sqrt((V2resp_centroid{2}(:,1).*50).^2+(V2resp_centroid{2}(:,2).*200).^2);
centshiftV2=sqrt((V2resp_centroid{6}(:,1).*50).^2+(V2resp_centroid{6}(:,2).*200).^2);

%roundedV1=round(ElectdistV1,-2);%nearest 100
roundedV1=round(ElectdistV1./2,-2).*2;%nearest 200

avFR=zeros(4,1);
stdFR=zeros(4,1);
avgnumelect=zeros(4,1);
stdavgnumelect=zeros(4,1);
avgsinglenumelect=zeros(4,1);
avFRsingle=zeros(4,1);
stdavgnumelectsingle=zeros(4,1);
stdFRsingle=zeros(4,1);
avgcent=zeros(4,1);
stdcent=zeros(4,1);
datforsig=nan(4,1000);
for i=200:200:800
    avgcent(i/200)=mean(centshiftV2(roundedV1==i),'omitnan');
    stdcent(i/200)=SEM(centshiftV2(roundedV1==i),0);
    avFR(i/200)=mean(V2resp_centroid{3}(roundedV1==i),'omitnan');
     avFRsingle(i/200)=mean(V2resp_centroid{5}(roundedV1==i),'omitnan');
    datforsig(i/200,1:sum(roundedV1==i))=V2resp_centroid{1}(roundedV1==i);
    stdFR(i/200)=SEM(V2resp_centroid{3}(roundedV1==i),0);
    stdFRsingle(i/200)=SEM(V2resp_centroid{5}(roundedV1==i),0);
        avgnumelect(i/200)=mean(V2resp_centroid{1}(roundedV1==i),'omitnan');
        avgsinglenumelect(i/200)=mean(V2resp_centroid{4}(roundedV1==i),'omitnan');
    stdavgnumelect(i/200)=SEM(V2resp_centroid{1}(roundedV1==i),0);
       stdavgnumelectsingle(i/200)=SEM(V2resp_centroid{4}(roundedV1==i),0);
    sum(~isnan(V2resp_centroid{5}(roundedV1==i)))
end

figure
hold on
%errorbar(100:100:800,avgcent,stdcent)
errorbar(200:200:800,avFR.*1000,stdFR.*1000)
xlabel('V1 electrode separation distance(\mum)')
ylabel('Firing rate (sp/s)')
beautifyPlot;
set(gca,'TickDir','out');
title('RFs')
ylim([0.4 1.3])
xlim([100 800])

figure
hold on
%errorbar(100:100:800,avgcent,stdcent)
errorbar(200:200:800,avgnumelect,stdavgnumelect)
xlabel('V1 electrode separation distance(\mum)')
ylabel('Number active electrodes in V2')
beautifyPlot;
set(gca,'TickDir','out');
title('RFs')
xlim([100 800])

figure
hold on
%errorbar(100:100:800,avgcent,stdcent)
hold on
errorbar(200:200:800,avFR.*1000,stdFR.*1000)
errorbar(200:200:800,avFRsingle.*1000,stdFRsingle.*1000)
xlabel('V1 electrode separation distance(\mum)')
ylabel('Firing rate (sp/s)')
beautifyPlot;
set(gca,'TickDir','out');
title('RFs')
xlim([100 800])

figure
hold on
%errorbar(100:100:800,avgcent,stdcent)
errorbar(200:200:800,avgnumelect,stdavgnumelect)
errorbar(200:200:800,avgsinglenumelect,stdavgnumelectsingle)
xlabel('V1 electrode separation distance(\mum)')
ylabel('Number active electrodes in V2')
beautifyPlot;
set(gca,'TickDir','out');
title('RFs')
xlim([100 800])

figure
hold on
%errorbar(100:100:800,avgcent,stdcent)
errorbar(200:200:800,avgcent,stdcent)

xlabel('V1 electrode separation distance(\mum)')
ylabel('Centroid shift (\mum)')
beautifyPlot;
set(gca,'TickDir','out');
title('RFs')
xlim([100 800])

%% Centroid
% timetolook=92:181;
% s = RandStream('mt19937ar','Seed',296);
% subselectbaseline=1:89;
% permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
% subselectbaseline=permuted_numbers(1:50);
% 
% E_MAP=Depth(1);
% emap_shaped=reshape(E_MAP,16,8);
% emap_shaped=emap_shaped(:,[1 3 4 2 5 7 8 6]);
% 
% numelecttotal=length(savefilename{1}{2}.CHN);
% V2resp_centroid=cell(numelecttotal,1);
% V2_centroid_array=cell(numelecttotal,1);
% for numstimelect=1:numelecttotal
%     it=0;
%     V2_centroid_array{numstimelect}=nan(31,7,300);
%     for numerfolders=1:numfolderstotal
%         trialend=size(savefilename{numerfolders}{2}.AMP,1);
%         for trial=1:trialend
%             whichelectstim=savefilename{numerfolders}{2}.AMP(trial,:)==ampinterest;
%             if sum(whichelectstim)~=numstimelect%check if the trial has the correct num stim elect and correct amplitude stim
%                 continue
%             end
%           
%             
%             V2only=mean(ratespiking{numerfolders}(1:64,timetolook,trial)-mean(ratespiking{numerfolders}(1:64,subselectbaseline,trial),2,'omitnan'),2,'omitnan');
%             if any(V2only>0)
%             [centV2r,centV2c]=centroidpos(V2only(emap_shaped(:,1:4)));
%             V2resp_centroid{numstimelect}=cat(1,V2resp_centroid{numstimelect},[centV2r,centV2c]);
%             it=it+1;
%             V2_centroid_array{numstimelect}(16-round(centV2r)+1:16-round(centV2r)+16,4-round(centV2c)+1:4-round(centV2c)+4,it)=V2only(emap_shaped(:,1:4)).*1000;
%             end
%         end
%     end
%     figure(numstimelect+100)
%     hold on
%     imagesc(mean(V2_centroid_array{numstimelect},3,'omitnan'))
%     title(num2str(numstimelect))
%     
% end



%% rate psth
% figure
% ax=axes;
% hold on
% color1 = linspace(0,1,numelecttotal);
% newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
% legendCell=cell(numelecttotal,1);
% for numelect=1:numelecttotal
%    stdshade(V1resp{numelect}.*1000,0.2,newcolors(numelect,:),[-90:90],1,ax); 
%    legendCell{numelect} = num2str(numelect);
% end
% xlim([-50 85])
% leg=legend(legendCell);
% title(leg,'# stim elect')
% xlabel('Time (ms)')
% ylabel('Firing rate (sp/s)')
% title('V1')
% beautifyPlot;
% 
% 
% figure
% ax=axes;
% hold on
% color1 = linspace(0,1,numelecttotal);
% newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
% legendCell=cell(numelecttotal,1);
% for numelect=1:numelecttotal
%    stdshade(V2resp{numelect}.*1000,0.2,newcolors(numelect,:),[-90:90],1,ax); 
%    legendCell{numelect} = num2str(numelect);
% end
% xlim([-50 85])
% leg=legend(legendCell);
% title(leg,'# stim elect')
% xlabel('Time (ms)')
% ylabel('Firing rate (sp/s)')
% title('V2')
% beautifyPlot;
% 
% %%
% for whichelectfocus=1:numelecttotal
% figure
% ax=axes;
% hold on
% color1 = linspace(0,1,numelecttotal);
% newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
% legendCell=cell(numelecttotal,1);
% for numelect=1:numelecttotal
%    stdshade(separateeachstimelect{whichelectfocus,numelect}.*1000,0.2,newcolors(numelect,:),[-90:90],1,ax); 
%    legendCell{numelect} = num2str(numelect);
% end
% xlim([-50 85])
% leg=legend(legendCell);
% title(leg,'# stim elect')
% xlabel('Time (ms)') 
% ylabel('Firing rate (sp/s)')
% title(['V2 elect ' num2str(whichelectfocus)])
% beautifyPlot;
% end
%%
% 2. # elect vs spread

% 3. # elect vs distribution early late, large small



end