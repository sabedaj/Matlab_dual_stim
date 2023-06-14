function epochratespike(ratestruct,savefilename,chnrange,excitesupress, chnexclude,normalisedat)
%% epoch has to be increasing with current once significant
% trial avg single elect - - standard deviation -  bins all electrodes based on thresh
if normalisedat==1
    multiplyspk=1;
else
    multiplyspk=1000;
end
D_data=dir;
if length(D_data)>size(ratestruct)
cd([D_data(12).folder filesep D_data(12).name;])
end
ratio=nan(1000,1);
E_MAP=Depth(1);
ordershapearray=reshape(E_MAP,16,8);
ordershapearray=ordershapearray(:,[1,3,4,2]);
alldata=cell(length(savefilename{15}{2}.AMP),1);
trialsig=cell(9,length(savefilename{15}{2}.AMP),length(savefilename));
folderIDsig=cell(9,length(savefilename{15}{2}.AMP));
stimchncount=0;
spread=nan(length(savefilename{15}{2}.AMP),500);
spreadgroup=nan(length(savefilename{15}{2}.AMP),9,500);
groupdata=cell(length(savefilename{15}{2}.AMP),1);
correlationdata=cell(length(savefilename{15}{2}.AMP),1);
stimchnsignificant=cell(length(savefilename{15}{2}.AMP),1);
allV1electsignificantstimchn=cell(length(savefilename{15}{2}.AMP),1);
alldata10uA=cell(5,1);
 totalchncount=zeros(5,1);
heatmap_centroid=cell(5,1);
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=cell(numfolderstotal,1);
for ampit=1:length(savefilename{15}{2}.AMP)
    heatmap_centroid{ampit}=nan(31,7,500);
ampinterest=savefilename{15}{2}.AMP(ampit);
groupdata{ampit}=cell(9,1);
correlationdata{ampit}=cell(9,1);
stimchnsignificant{ampit}=cell(9,1);
allV1electsignificantstimchn{ampit}=cell(9,1);
alldata10uA{ampit}=cell(9,1);
significantspread=nan(128,500);
stimchncount=0;
significantspread_groupsplit=nan(128,500,9);

for numerfolders=1:numfolderstotal

    trialend=length(fieldnames(ratestruct{numerfolders}));
    for trial=1:trialend
        trialcheck=['T' num2str(trial)];
        ratespiking{numerfolders}(:,:,trial)=nanmean(ratestruct{numerfolders}.(trialcheck),3);
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
   rateAMPua(chnexclude,:,:)=nan;
   if normalisedat==1
   rateAMPua=rateAMPua./max(squeeze(max(rateAMPua,[],2)),[],2);
   end
   simchnall=savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),3);
   stimchnall_notsorted=E_MAP(simchnall);
   totalchncount(ampit)=length(chnrange)*size(rateAMPua,3)+totalchncount(ampit);
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
   %edit savefilename{4} in loop to store channels in 3rd column. then can
   %use the trial IDs from
   %savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1)
   %to call individual channels and store in stimchnsignificant -> make sure
   %only stored once though
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %stimchnall=;
CatergoriesEBL=cell(size(sig,3),1);
   for chn=chnrange
       for stimchn=1:size(sig,3)
           if isempty(CatergoriesEBL{stimchn})
               CatergoriesEBL{stimchn}=cell(3,1);
           end
           if check_multiple(chn,stimchn)
               group=find(squeeze(sig(chn,:,stimchn)==1));
               if ~isempty(group) && all(group<5) 
                   CatergoriesEBL{stimchn}{1}=[CatergoriesEBL{stimchn}{1}; rateAMPua(chn,:,stimchn)];
               elseif ~isempty(group) && all(group>5) %&& all(group<9)
                   CatergoriesEBL{stimchn}{3}=[CatergoriesEBL{stimchn}{3}; rateAMPua(chn,:,stimchn)];
               elseif ~isempty(group) && any(group>5) && any(group<5)% && all(group<9)
                   CatergoriesEBL{stimchn}{2}=[CatergoriesEBL{stimchn}{2}; rateAMPua(chn,:,stimchn)];
               end
               sig(chn,:,stimchn)=false;
               sig(chn,pmax(chn,stimchn),stimchn)=true;
               group=find(squeeze(sig(chn,:,stimchn)==1));
           else
               group=find(squeeze(sig(chn,:,stimchn)==1));
               if ~isempty(group) && group<5 
                   CatergoriesEBL{stimchn}{1}=[CatergoriesEBL{stimchn}{1}; rateAMPua(chn,:,stimchn)];
               elseif ~isempty(group) && all(group>5) %&& all(group<9)
                   CatergoriesEBL{stimchn}{3}=[CatergoriesEBL{stimchn}{3}; rateAMPua(chn,:,stimchn)];
               end
           end
          
           if ~isempty(group)
               for gnum=1:length(group)
                   groupdata{ampit}{group(gnum)}=[groupdata{ampit}{group(gnum)}; rateAMPua(chn,:,stimchn)];
 
                   alldata{ampit}=[alldata{ampit}; rateAMPua(chn,:,stimchn)];
                   [co,~]=xcorr(rateAMPua(stimchnall_notsorted(stimchn),90:180,stimchn).*multiplyspk,rateAMPua(chn,90:180,stimchn).*multiplyspk,'coeff');
                   correlationdata{ampit}{group(gnum)}=[correlationdata{ampit}{group(gnum)}; co];
                   
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   if (isempty(stimchnsignificant{ampit}{group(gnum)}) || ~any(all(stimchnsignificant{ampit}{group(gnum)}==rateAMPua(stimchnall_notsorted(stimchn),:,stimchn),2)))
                       stimchnsignificant{ampit}{group(gnum)}=[stimchnsignificant{ampit}{group(gnum)}; rateAMPua(stimchnall_notsorted(stimchn),:,stimchn)];
                       if stimchnall_notsorted(1)>64
                           sigV1=squeeze(any(sig,2));
                           sigV1(1:64,:)=false;
                           allV1electsignificantstimchn{ampit}{group(gnum)}=[allV1electsignificantstimchn{ampit}{group(gnum)}; rateAMPua(sigV1(:,stimchn),:,stimchn)];
                       else
                           allV1electsignificantstimchn{ampit}{group(gnum)}=[allV1electsignificantstimchn{ampit}{group(gnum)}; rateAMPua(1:64,:,stimchn)];
                       end
                   end
                   %use if you only want the channels in 10uA from 2:10uA
                   if ampit==5
                       for itAua=1:5
                           ampinterest1=savefilename{15}{2}.AMP(itAua);
                           [m,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest1)); ampinterest1=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp
                           
                           rateAMPua1=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest1),1));
                           alldata10uA{itAua}{group(gnum)}=[alldata10uA{itAua}{group(gnum)}; rateAMPua1(chn,:,stimchn)];
                       end
                   end
                   significantspread(chn,stimchn+stimchncount)=meancatdata(chn,group(gnum),stimchn);
                   significantspread_groupsplit(chn,stimchn+stimchncount,group(gnum))=meancatdata(chn,group(gnum),stimchn);
                   trialsig{group(gnum),ampit,numerfolders}=[trialsig{group(gnum),ampit,numerfolders} stimchn];%determine how many different stim electrodes creating significance
                   folderIDsig{group(gnum),ampit}=[folderIDsig{group(gnum),ampit} str2double(savefilename{numerfolders}{3})];%determine how many different monkeys creating significance
               end
           end
           
       end
       
   end
    
   figure(20*ampit)
   hold on
   for i=1:size(sig,3)
       if (size(CatergoriesEBL{i}{1},1)+size(CatergoriesEBL{i}{3},1))>2
       ratio(stimchncount+i)=size(CatergoriesEBL{i}{3},1)/(size(CatergoriesEBL{i}{1},1)+size(CatergoriesEBL{i}{3},1));
        scatter(stimchncount+i,ratio(stimchncount+i),'k')
       end
   end
   ylabel('Early(0) vs late(1)')
   xlabel('Stimchn session #')
  stimchncount=stimchncount+size(sig,3);
    
end
fprintf(num2str(stimchncount))
%current significance
orderedsigchns=significantspread(E_MAP,:);
avgspread=nan(size(orderedsigchns,2),1);
figure(21+ampit)
hold on
hist(ratio)
   xlabel('Early(0) vs late(1)')
   ylabel('# stimchn sessions')
   title([num2str(savefilename{15}{2}.AMP(ampit)),'\muA'])
for itstimchn=1:size(orderedsigchns,2)
    tmp=reshape(orderedsigchns(chnrange,itstimchn),16,4);
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
        tmp=reshape(squeeze(orderedgroupsig(chnrange,itstimchn,group)),16,4);
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

%this merges all groups
orderedgroupsig=sum(significantspread_groupsplit(E_MAP,:,:),3,'omitnan');%merges the groups because we don't care about that for now
for itstimchn=1:size(orderedsigchns,2)
    tmp=reshape(orderedgroupsig(chnrange,itstimchn),16,4);
    tmp=tmp(:,[1 3 4 2]);
    tmp(isnan(tmp))=0;
    xcent=sum(tmp.*[1 2 3 4],'all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location
    ycent=sum(tmp.*(1:16)','all','omitnan')./sum((tmp),'all','omitnan');%centroid/centre of mass location

    if ~isnan(xcent) && ~isnan(ycent)
        xcent=round(xcent);
         ycent=round(ycent);
        heatmap_centroid{ampit}(16-ycent+1:16-ycent+16,4-xcent+1:4-xcent+4,itstimchn)=tmp;
    end
end

spreadgroup(ampit,1:9,1:length(avgspread))=avgspreadwg;
spread(ampit,1:length(avgspread))=avgspread;
figure(1000*ampit); hold on;
tmpdat=cellfun(@(x) mean(x,1,'omitnan').*multiplyspk, groupdata{ampit}, 'UniformOutput', false);
%use to remove lines of data
%tmpdat(2:2:end)=[];
cellfun(@(x) plot(-90:90,x), tmpdat)
color1 = linspace(0,1,size(tmpdat,1));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
indvsig=cellfun(@(x) size(x,1), groupdata{ampit});
text(-80,24,'# sig:')
text(-80,20,num2str(indvsig))
%text(-80,10,num2str(indvsig))V2
totalsig=sum(indvsig);
text(-80,40,['Total sig: ' num2str(totalsig) ' / ' num2str(totalchncount(ampit))])
xlabel('Time (ms)')
ylabel('Firing rate (Sp/s)')
xlim([-85 85])
if normalisedat==1
    ylim([0 1])
else
ylim([0 20])
end
%ylim([0 400])
set(gca,'TickDir','out');
title([num2str(ampinterest) '\muA'])
leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','62:71','72:81','82:91','Average');
title(leg,'Time epoch(ms)')
figure (6)
plot(indvsig)
hold on
xlabel('V2 Epochs')
ylabel('# channels significantly responding in V2')



figure(1001*ampit); hold on;tmpdat=cellfun(@(x) mean(x,1,'omitnan').*multiplyspk, stimchnsignificant{ampit}, 'UniformOutput', false);
%use to remove lines of data
%tmpdat(2:2:end)=[];
cellfun(@(x) plot(-90:90,x), tmpdat)
color1 = linspace(0,1,size(tmpdat,1));
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);

indvsig=cellfun(@(x) size(x,1), stimchnsignificant{ampit});
text(-80,350,['# sig / ' num2str(stimchncount) ' :'])
text(-80,250,num2str(indvsig))
xlabel('Time (ms)')
ylabel('Firing rate (Sp/s)')
xlim([-85 85])
%ylim([0 25])
if normalisedat==1
    ylim([0 1])
else
ylim([0 200])
end
set(gca,'TickDir','out');
title([num2str(ampinterest) '\muA Stimchn'])



figure(1002*ampit); hold on;cellfun(@(x) plot(-90:90,mean(x,1,'omitnan').*multiplyspk), allV1electsignificantstimchn{ampit})
color1 = linspace(0,1,groupit);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
indvsig=cellfun(@(x) size(x,1), allV1electsignificantstimchn{ampit});
text(-80,350,['# sig / ' num2str(stimchncount*(64-15)) ' :'])
text(-80,250,num2str(indvsig))
xlabel('Time (ms)')
ylabel('Firing rate (Sp/s)')
xlim([-85 85])
%ylim([0 25])
ylim([0 400])
set(gca,'TickDir','out');
title([num2str(ampinterest) '\muA All chns on stim array'])
figure (5)
plot(indvsig)
hold on
xlabel('V2 Epochs')
ylabel('# channels significantly responding in V1')

figure(1003*ampit); hold on;cellfun(@(x) plot(-90:90,mean(x,1,'omitnan')), correlationdata{ampit})
color1 = linspace(0,1,groupit);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
xlabel('Lag (ms)')
ylabel('Correlation(Sp/s)')
set(gca,'TickDir','out');
title(['Correlation' num2str(ampinterest) '\muA'])
leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','62:71','72:81','82:91','Average');
title(leg,'Time epoch(ms)')

% if ampit==5 % plots V1 and V2 together for each epoch individually 
%     for epoch=1:9
% figure
% title(['Epoch ' num2str(epoch) ' ' num2str(ampinterest) '\muA'])
% hold on
% plot(-90:90,mean(stimchnsignificant{ampit}{epoch},'omitnan').*multiplyspk)
% plot(-90:90,mean(groupdata{ampit}{epoch},'omitnan').*multiplyspk)
% xlabel('Time (ms)')
% ylabel('Firing rate (Sp/s)')
% xlim([-85 85])
% if normalisedat==1
%     ylim([0 1])
% else
% ylim([0 200])
% end
% set(gca,'TickDir','out');
%     end
% 
% end
% figure(15);hold on
% plot(1:3,cellfun(@(x) size(x,1),CatergoriesEBL));
% xlabel('Categories')
% xticks([1 2 3])
% xticklabels(gca, {'Early','Both','Late'});
% set(gca,'TickDir','out');
% ylabel('# V2 electrodes')
end
%% peak and time to return to baseline
%V1
maxpeakV1=cellfun(@(x) max(x(:,90:end),[],2),stimchnsignificant{5},'UniformOutput',false);
avgpeakV1=cellfun(@(x) mean(x,'omitnan').*multiplyspk,maxpeakV1);

%V2
maxpeakV2=cellfun(@(x) max(x(:,90:end),[],2),groupdata{5},'UniformOutput',false);
avgpeakV2=cellfun(@(x) mean(x,'omitnan').*multiplyspk,maxpeakV2);
figure; plot(1:9,avgpeakV1)
hold on
plot(1:9,avgpeakV2)
set(gca,'TickDir','out');
xlabel('epoch')
ylabel('Max firing rate')

leg=legend('V1','V2')

%width resp
timemaxsave=nan(300,9);
for group=1:9
    for i = 1:size(maxpeakV1{group})
        if ~isnan(maxpeakV1{group}(i)) && maxpeakV1{group}(i)~=0
            [~,c]=find(stimchnsignificant{5}{group}(i,90:end)./maxpeakV1{group}(i)==1,1,'first');
            baselinemean=mean(stimchnsignificant{5}{group}(i,1:85),'omitnan');
            baselinestd=std(stimchnsignificant{5}{group}(i,1:85),'omitnan');
            [~,c2]=find(stimchnsignificant{5}{group}(i,c:end)<=baselinemean+baselinestd,1,'first');
            timemaxsave(i,group)=c+c2;
        end
    end
end
timemax=timemaxsave;
figure; plot(1:9, mean(timemax,'omitnan'))
set(gca,'TickDir','out');
xlabel('epoch')
ylabel('Latency to return to baseline after peak (ms)')
title('V1')


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


figure; heatmap(sum(heatmap_centroid{5}>0,3)); set(gca,'ColorScaling','log')
dataheatmap=sum(heatmap_centroid{5}>0,3);
figure; heatmap(dataheatmap,'CellLabelColor','none','GridVisible','off');set(gca,'ColorScaling','log')
maxdata=max(dataheatmap,[],'all')/2;
fwhmdat=dataheatmap>=maxdata;



%% sigmoid using groupdata
sigmoidampepoch=zeros(5,9);
for ampit=1:5
    for epoch=1:9
        sigmoidampepoch(ampit,epoch)=mean(alldata10uA{ampit}{epoch}(:,82+(10*epoch):91+(10*epoch)),'all'); 
    end
end
figure; hold on
for i=1:9
plot([2 5 6 8 10],sigmoidampepoch(:,i).*multiplyspk)
end
color1 = linspace(0,1,9);
newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
colororder(newcolors);
%plot([2 5 6 8 10],mean(sigmoidampepoch,2).*1000,'k')
xlabel('current (\muA)')
ylabel('Firing rate (sp/s)')
leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','61:72','72:81','82:91','Average');
%leg=legend('2:11','22:31','42:51','61:72','82:91');
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
%leg=legend('2:11','22:31','42:51','62:71','82:91','Total','Significance');
leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','61:72','72:81','82:91','Total','Significance');
title(leg,'Time epoch(ms)')
set(gca,'TickDir','out');

end
