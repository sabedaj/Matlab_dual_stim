function epochratespike(ratestruct,savefilename,chnrange,excitesupress, chnexclude,normalisedat)
%% epoch has to be increasing with current once significant
% trial avg single elect - - standard deviation -  bins all electrodes based on thresh
if normalisedat==1
    multiplyspk=1;
else
    multiplyspk=1000;
end
 s = RandStream('mt19937ar','Seed',296);%45
D_data=dir;
if length(D_data)>size(ratestruct)
cd([D_data(12).folder filesep D_data(12).name;])
end
ratio=nan(1000,1);
E_MAP=Depth(1);
centredstimchn=cell(5,1);
stimchnsave=cell(5,1);
foldercheck=cell(5,1);
ordershapearray=reshape(E_MAP,16,8);
ordershapearray=ordershapearray(:,[1,3,4,2]);
alldata=cell(length(savefilename{15}{2}.AMP),1);
alldata_stim=cell(length(savefilename{15}{2}.AMP),1);

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
V2ELresp=cell(5,1);
V1ELresp=cell(5,1);
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=cell(numfolderstotal,1);
for ampit=1:length(savefilename{15}{2}.AMP)
    iteratestimdat=1;
iteratealldat=1;
    alldata_stim{ampit}=nan(108000,181);
    alldata{ampit}=nan(1800,181);
    centredstimchn{ampit}=nan(31,7,500);
    heatmap_centroid{ampit}=nan(31,7,500);
ampinterest=savefilename{15}{2}.AMP(ampit);
V2ELresp{ampit}=cell(3,1);
V1ELresp{ampit}=cell(3,1);
groupdata{ampit}=cell(9,1);
correlationdata{ampit}=cell(9,1);
stimchnsignificant{ampit}=cell(9,1);
allV1electsignificantstimchn{ampit}=cell(9,1);
alldata10uA{ampit}=cell(9,1);
significantspread=nan(128,500);
stimchncount=0;
significantspread_groupsplit=nan(128,500,9);
foldercheck{ampit}=cell(numfolderstotal,1);
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
       rateELcurrent=nan(128,2,length(savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1)),5);
       rateecurrent=nan(128,length(savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1)),5);

       for ampit1=1:length(savefilename{15}{2}.AMP)
           ampinterest1=savefilename{15}{2}.AMP(ampit1);
           [m,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest1)); ampinterest1=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp
           
           rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest1),1));
           rateepochcurrent(1:size(rateAMPua,1),:,:,ampit1)=cat(2,mean(rateAMPua(:,92:101,:),2),mean(rateAMPua(:,102:111,:),2),mean(rateAMPua(:,112:121,:),2),mean(rateAMPua(:,122:131,:),2),mean(rateAMPua(:,132:141,:),2),mean(rateAMPua(:,142:151,:),2),mean(rateAMPua(:,152:161,:),2),mean(rateAMPua(:,162:171,:),2),mean(rateAMPua(:,172:181,:),2));
           rateELcurrent(1:size(rateAMPua,1),:,:,ampit1)=cat(2,mean(rateAMPua(:,92:136,:),2),mean(rateAMPua(:,137:181,:),2));
           rateecurrent(1:size(rateAMPua,1),:,ampit1)=mean(rateAMPua(:,92:181,:),2);

       end
       validchannels=sum(diff(rateepochcurrent,[],4),4,'omitnan')>0;
       validchannelsEL=sum(diff(rateELcurrent,[],4),4,'omitnan')>0;
       validchannelsingle=sum(diff(rateecurrent,[],3),3,'omitnan')>0;

   end
   rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1));
   rateAMPua(chnexclude,:,:)=nan;
   if normalisedat==1
   rateAMPua=rateAMPua./max(squeeze(max(rateAMPua,[],2)),[],2);
   end

   goodchns=chnrange(~ismember(chnrange,unique(chnexclude)));
   totalchncount(ampit)=(length(goodchns))*size(rateAMPua,3)+totalchncount(ampit);
   meancatdata=cat(2,mean(rateAMPua(:,92:101,:),2),mean(rateAMPua(:,102:111,:),2),mean(rateAMPua(:,112:121,:),2),mean(rateAMPua(:,122:131,:),2),mean(rateAMPua(:,132:141,:),2),mean(rateAMPua(:,142:151,:),2),mean(rateAMPua(:,152:161,:),2),mean(rateAMPua(:,162:171,:),2),mean(rateAMPua(:,172:181,:),2));
   meancatdataEL=cat(2,mean(rateAMPua(:,92:136,:),2),mean(rateAMPua(:,137:181,:),2));
   sig=false(128,9,size(meancatdata,3));
   sigEL=false(128,2,size(meancatdata,3));
   if excitesupress==0
       thresh=squeeze(mean(rateAMPua(:,1:85,:),2)-(std(rateAMPua(:,1:85,:),[],2).*1)); %%%%%NOTE LOWER THRESHOLD
       for groupit=1:9
           sig(:,groupit,:)=squeeze(meancatdata(:,groupit,:))<thresh;
       end
       check_multiple=sum(sig,2)>1;
       [~,pmax]=min(meancatdata,[],2);%pmin
       pmax=squeeze(pmax);%min
   else
       subselectbaseline=1:89;
       permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
       subselectbaseline=permuted_numbers(1:50);
       thresh=squeeze(mean(rateAMPua(:,subselectbaseline,:),2)+(std(rateAMPua(:,subselectbaseline,:),[],2).*3));
       for groupit=1:9
           sig(1:size(meancatdata,1),groupit,:)=squeeze(meancatdata(:,groupit,:))>thresh;
       end
       sig=sig&validchannels;
       %permutedRate = permute(rateAMPua, [2 1 3]);
       sigsingle=squeeze(mean(rateAMPua(:,92:181,:),2))>thresh;%ttest(permutedRate(92:181,:,:),permutedRate(1:90,:,:),"Tail","right");%
       %sigsingle(isnan(sigsingle))=false;
       %sigsingle=squeeze(sigsingle)&validchannelsingle;
       check_multiple=squeeze(sum(sig,2)>1);
       [~,pmax]=max(meancatdata,[],2);
       pmax=squeeze(pmax);
       for EL=1:2
            sigEL(1:size(meancatdataEL,1),EL,:)=squeeze(meancatdataEL(:,EL,:))>thresh;
       end
       sigEL=sigEL&validchannelsEL;
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
            groupEL=find(squeeze(sigEL(chn,:,stimchn)==1));
           if ~isempty(groupEL) && length(groupEL)==2
               CatergoriesEBL{stimchn}{2}=[CatergoriesEBL{stimchn}{2}; rateAMPua(chn,:,stimchn)];%both
           elseif ~isempty(groupEL) && groupEL<2
               CatergoriesEBL{stimchn}{1}=[CatergoriesEBL{stimchn}{1}; rateAMPua(chn,:,stimchn)];%early
           elseif ~isempty(groupEL) && all(groupEL>1) %&& all(group<9)
               CatergoriesEBL{stimchn}{3}=[CatergoriesEBL{stimchn}{3}; rateAMPua(chn,:,stimchn)];%late
           end
           
           if check_multiple(chn,stimchn)
               sig(chn,:,stimchn)=false;
               sig(chn,pmax(chn,stimchn),stimchn)=true;
           end
           simchnall=savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),3);
           stimchnall_notsorted=E_MAP(simchnall);
           if sigsingle(chn,stimchn)==1
               alldata{ampit}(iteratealldat,1:181)=rateAMPua(chn,:,stimchn);
               iteratealldat=iteratealldat+1;
               %if  ~any(foldercheck{ampit}{numerfolders}==stimchn)
               alldata_stim{ampit}(iteratestimdat:iteratestimdat+63,1:181)=rateAMPua(65:128,:,stimchn);
               iteratestimdat=iteratestimdat+64;
               significantspread(chn,stimchn+stimchncount)=mean(rateAMPua(chn,92:181,stimchn),2);
               folderIDsig{ampit}=[folderIDsig{ampit} str2double(savefilename{numerfolders}{3})];%determine how many different monkeys creating significance
               %stimchnsave{ampit}=[stimchnsave{ampit}; stimchn];
               %foldercheck{ampit}{numerfolders}=[foldercheck{ampit}{numerfolders};stimchn];
               % end
              
                if all(chnrange>64)
                    stimchnarray=[1:16;33:48;49:64;17:32]';
          
                     sortedmap=reshape(E_MAP(65:128),16,4);
                     sortedmap=sortedmap(:,[1 3 4 2]);
                     [rchn,cchn]=find(E_MAP(chn)==sortedmap);
                    %V1array{ampit}(chn,:)=squeeze(mean(rateAMPua(chn,92:181,:),2));
                    
                     [r,c]=find(stimchnarray==(simchnall(stimchn)-64));
                     centredstimchn{ampit}(17-r+(rchn-1),5-c+(cchn-1),stimchncount+stimchn)=squeeze(mean(rateAMPua(chn,92:102,stimchn),2));
%                     V1stimelectposRC{ampit}(stimchn,:)=simchnall
%                     V1arraysorted=reshape(V1array(E_MAP(1:64),:),[16,4,size(V1array,2)]);
%                     V1arraysorted=V1arraysorted(:,[1 3 4 2],:);
%                     for i=1:length(simchnall)
%                         [r,c]=find(stimchnarray==(simchnall(i)-64));
%                         centredstimchn{ampit}(17-r:32-r,5-c:8-c,stimchncount+i)=V1arraysorted(:,:,i);
%                     end
                    
                end
          end
          group=find(squeeze(sig(chn,:,stimchn)==1));
           if ~isempty(group)
               for gnum=1:length(group)
                   groupdata{ampit}{group(gnum)}=[groupdata{ampit}{group(gnum)}; rateAMPua(chn,:,stimchn)];
 
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
                 
                   significantspread_groupsplit(chn,stimchn+stimchncount,group(gnum))=meancatdata(chn,group(gnum),stimchn);
                   trialsig{group(gnum),ampit,numerfolders}=[trialsig{group(gnum),ampit,numerfolders} stimchn];%determine how many different stim electrodes creating significance
                   %folderIDsig{group(gnum),ampit}=[folderIDsig{group(gnum),ampit} str2double(savefilename{numerfolders}{3})];%determine how many different monkeys creating significance
               end
           end
           
       end
       
   end
    
%    figure(20*ampit)
%    hold on
   for i=1:size(sig,3)
       if (size(CatergoriesEBL{i}{1},1)+size(CatergoriesEBL{i}{3},1))>0
       ratio(stimchncount+i)=size(CatergoriesEBL{i}{3},1)/(size(CatergoriesEBL{i}{1},1)+size(CatergoriesEBL{i}{3},1));%determines the ratio of early to late responses in V2
      
       %scatter(stimchncount+i,ratio(stimchncount+i),'k')
       end
       %determine if the early and late responses have different
       %characteristics
        numcategories=cellfun(@(x) ~isempty(x), CatergoriesEBL{i});
        if sum(numcategories)==1
            V2ELresp{ampit}{numcategories}=[V2ELresp{ampit}{numcategories}; CatergoriesEBL{i}{numcategories}];
            V1ELresp{ampit}{numcategories}=[V1ELresp{ampit}{numcategories}; rateAMPua(65:128,:,stimchn)];
           
        end
   end
%    ylabel('Early(0) vs late(1)')
%    xlabel('Stimchn session #')
  stimchncount=stimchncount+size(sig,3);
    
end
fprintf(num2str(stimchncount))
%current significance
orderedsigchns=significantspread(E_MAP,:);
avgspread=nan(size(orderedsigchns,2),1);
% figure(21+ampit)
% hold on
% hist(ratio)
%    xlabel('Early(0) vs late(1)')
%    ylabel('# stimchn sessions')
%    title([num2str(savefilename{15}{2}.AMP(ampit)),'\muA'])
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
orderedgroupsig=significantspread(E_MAP,:);%merges the groups because we don't care about that for now
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
alldata{ampit}(iteratealldat:end,:)=[];
alldata_stim{ampit}(iteratestimdat:end,:)=[];
% spreadgroup(ampit,1:9,1:length(avgspread))=avgspreadwg;
% spread(ampit,1:length(avgspread))=avgspread;
% figure(1000*ampit); hold on;
% tmpdat=cellfun(@(x) mean(x,1,'omitnan').*multiplyspk, groupdata{ampit}, 'UniformOutput', false);
% %use to remove lines of data
% %tmpdat(2:2:end)=[];
% cellfun(@(x) plot(-90:90,x), tmpdat)
% color1 = linspace(0,1,size(tmpdat,1));
% newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
% colororder(newcolors);
% indvsig=cellfun(@(x) size(x,1), groupdata{ampit});
% text(-80,24,'# sig:')
% text(-80,20,num2str(indvsig))
% %text(-80,10,num2str(indvsig))V2
% totalsig=sum(indvsig);
% text(-80,40,['Total sig: ' num2str(totalsig) ' / ' num2str(totalchncount(ampit))])
% xlabel('Time (ms)')
% ylabel('Firing rate (Sp/s)')
% xlim([-85 85])
% if normalisedat==1
%     ylim([0 1])
% else
% ylim([0 20])
% end
% %ylim([0 400])
% set(gca,'TickDir','out');
% title([num2str(ampinterest) '\muA'])
% leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','62:71','72:81','82:91','Average');
% title(leg,'Time epoch(ms)')
% figure (6)
% plot(indvsig)
% hold on
% xlabel('V2 Epochs')
% ylabel('# channels significantly responding in V2')



if any(ampit==[1 3 5])
    figure(55)

if ampit==1
    ax=axes;
    hold on
xlabel('Time epoch(ms)')
ylabel('Firing rate (sp/s)')
title('All significant channels averaged')
text(-80,0.0025,['# sig ' num2str(size(alldata{ampit},1)) ' / ' num2str(totalchncount(ampit)) ' electrodes'])
text(-80,0.0023,['Need ' num2str(totalchncount(ampit)*0.003) 'to be sig'])
xlim([-85 85])
set(gca,'TickDir','out');
end
color1 = linspace(0,1,5);
newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
colororder(newcolors);
[lineOut, fillOut] = stdshade(alldata{ampit}.*multiplyspk,0.1,newcolors(ampit,:),[-90:90],1,ax);
leg=legend('2','6','10');
title(leg,'Current(uA)')
xlim([-50 85])
%plot(-90:90,mean(alldata{ampit}.*multiplyspk))
end

% figure(51)
% hold on
% errorstd=std(alldata{ampit}(:,92:181).*multiplyspk,0,'all','omitnan');
% errorbar(savefilename{15}{2}.AMP(ampit),mean(alldata{ampit}(:,92:181).*multiplyspk,'all'),errorstd,'o')
% xlabel('Current(\muA)')
% ylabel('Firing rate (sp/s)')
% title('All significant channels averaged')
% 
% figure(1001*ampit); hold on;tmpdat=cellfun(@(x) mean(x,1,'omitnan').*multiplyspk, stimchnsignificant{ampit}, 'UniformOutput', false);
% %use to remove lines of data
% %tmpdat(2:2:end)=[];
% cellfun(@(x) plot(-90:90,x), tmpdat)
% color1 = linspace(0,1,size(tmpdat,1));
% newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
% colororder(newcolors);
% 
% indvsig=cellfun(@(x) size(x,1), stimchnsignificant{ampit});
% text(-80,350,['# sig / ' num2str(stimchncount) ' :'])
% text(-80,250,num2str(indvsig))
% xlabel('Time (ms)')
% ylabel('Firing rate (Sp/s)')
% xlim([-85 85])
% %ylim([0 25])
% if normalisedat==1
%     ylim([0 1])
% else
% ylim([0 200])
% end
% set(gca,'TickDir','out');
% title([num2str(ampinterest) '\muA Stimchn'])



% figure(1002*ampit); hold on;cellfun(@(x) plot(-90:90,mean(x,1,'omitnan').*multiplyspk), allV1electsignificantstimchn{ampit})
% color1 = linspace(0,1,groupit);
% newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
% colororder(newcolors);
% indvsig=cellfun(@(x) size(x,1), allV1electsignificantstimchn{ampit});
% text(-80,350,['# sig / ' num2str(stimchncount*(64-15)) ' :'])
% text(-80,250,num2str(indvsig))
% xlabel('Time (ms)')
% ylabel('Firing rate (Sp/s)')
% xlim([-85 85])
% %ylim([0 25])
% ylim([0 400])
% set(gca,'TickDir','out');
% title([num2str(ampinterest) '\muA All chns on stim array'])
% figure (5)
% plot(indvsig)
% hold on
% xlabel('V2 Epochs')
% ylabel('# channels significantly responding in V1')

% figure(1003*ampit); hold on;cellfun(@(x) plot(-90:90,mean(x,1,'omitnan')), correlationdata{ampit})
% color1 = linspace(0,1,groupit);
% newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
% colororder(newcolors);
% xlabel('Lag (ms)')
% ylabel('Correlation(Sp/s)')
% set(gca,'TickDir','out');
% title(['Correlation' num2str(ampinterest) '\muA'])
% leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','62:71','72:81','82:91','Average');
% title(leg,'Time epoch(ms)')

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

%% high and low V2 response vs V1 response
%need to sort this out. not what was expected at the moment. double
%dipping=cause? how to get around this?
avgtime=92:180;
baselinetime=subselectbaseline;%baseline to subtract only 50random samples
subselectbaseline=1:89;


HLB_data=cell(5,1);
baslineSV2=nan(5,1);
baselineSstdV2=baslineSV2;
avgresp2V2=baslineSV2;
stdresp2V2=baslineSV2;
Earlylate_data=cell(5,1);
baslineLV2=nan(5,1);
baselineLstdV2=baslineSV2;
avgrespLV2=baslineSV2;
stdrespLV2=baslineSV2;
% bigcorall=cell(5,1);
% smallcorall=cell(5,1);
DataLateV2=cell(5,1);
DataEarlyV2=cell(5,1);
DataLateV1=cell(5,1);
DataEarlyV1=cell(5,1);

for ampit=1:5
HLB_data{ampit}=cell(3,1);
avgalldata=mean(alldata{ampit}(:,avgtime),2).*multiplyspk;
[~, sorted_indices] = sort(avgalldata, 'ascend');
halflength=ceil(length(sorted_indices)/2);
sorted_indicesstim=[];

Earlylate_data{ampit}=cell(3,1);
avgalldata_Early=mean(alldata{ampit}(:,92:136),2).*multiplyspk;
avgalldata_Late=mean(alldata{ampit}(:,137:181),2).*multiplyspk;
sortLate=avgalldata_Early<avgalldata_Late;
DataLateV2{ampit}=alldata{ampit}(sortLate,:).*multiplyspk;
DataEarlyV2{ampit}=alldata{ampit}(~sortLate,:).*multiplyspk;

sortLateV1=repelem(sortLate,64);
DataLateV1{ampit}=alldata_stim{ampit}(sortLateV1,:).*multiplyspk;
DataEarlyV1{ampit}=alldata_stim{ampit}(~sortLateV1,:).*multiplyspk;
ttest(mean(DataLateV2{ampit}(:,92:181)-mean(DataLateV2{ampit}(:,baselinetime),2,'omitnan'),'omitnan'),mean(DataEarlyV2{ampit}(:,92:181)-mean(DataEarlyV2{ampit}(:,baselinetime),2,'omitnan'),'omitnan'))
ttest(mean(DataLateV1{ampit}(:,92:112)-mean(DataLateV1{ampit}(:,baselinetime),2,'omitnan'),'omitnan'),mean(DataEarlyV1{ampit}(:,92:112)-mean(DataEarlyV1{ampit}(:,baselinetime),2,'omitnan'),'omitnan'))

for i=1:length(sorted_indices)
sorted_indicesstim=[sorted_indicesstim sorted_indices(i)*64-63:sorted_indices(i)*64];
end
halflengthstim=halflength*64;

data=reshape(mean(alldata_stim{ampit}(sorted_indicesstim,:),2,'omitnan'),64,length(avgalldata));
data(isnan(data))=0;
[~,ia,~]=unique(data','rows');

alldatastim_sorted=alldata_stim{ampit}(sorted_indicesstim,:);
for iterateunique=1:length(ia)
    columncompare=data(:,ia(iterateunique));
    columnindex=find(ismember(data',columncompare','rows'));
    if all(columnindex>halflength)%big

        HLB_data{ampit}{1}=[HLB_data{ampit}{1}; alldatastim_sorted(columnindex(1)*64-63:columnindex(1)*64,:)];
    
      
%         for coriterate=1:length(columnindex)
%             stimdat=alldatastim_sorted(columnindex(coriterate)*64-63:columnindex(coriterate)*64,:);
%            
%                 for iteratestimarray=1:size(stimdat,1)
%                     [bigcor,~]=xcorr(alldata{ampit}(columnindex(coriterate),:)./max(alldata{ampit}(columnindex(coriterate),:)),stimdat(iteratestimarray,:)./max(stimdat(iteratestimarray,:)));
%                     bigcorall{ampit}=[bigcorall{ampit}; bigcor];
%                 end
%           
%         end
    elseif all(columnindex<halflength)%small
        HLB_data{ampit}{2}=[HLB_data{ampit}{2}; alldatastim_sorted(columnindex(1)*64-63:columnindex(1)*64,:)];
        %this kills 30mins to run - maybe try normalising by baseline? is
        %this valid??? what question are you asking> maybe normalise V2 but
        %not V1 - see if amplitude in v1 is correlated without influence fo
        %V2 mathematically. 
%         for coriterate=1:length(columnindex)
%             stimdat=alldatastim_sorted(columnindex(coriterate)*64-63:columnindex(coriterate)*64,:);
%             
%             for iteratestimarray=1:size(stimdat,1)
%                 [smallcor,~]=xcorr(alldata{ampit}(columnindex(coriterate),:)./max(alldata{ampit}(columnindex(coriterate),:)),stimdat(iteratestimarray,:)./max(stimdat(iteratestimarray,:)));
%                 smallcorall{ampit}=[smallcorall{ampit}; smallcor];
%             end
%             
%         end
    else%both
         HLB_data{ampit}{3}=[HLB_data{ampit}{3}; alldatastim_sorted(columnindex(1)*64-63:columnindex(1)*64,:)];
    end
end
HLB_data{ampit}{1}(sum(HLB_data{ampit}{1},2)==0,:)=[];
HLB_data{ampit}{2}(sum(HLB_data{ampit}{2},2)==0,:)=[];

%sigmoidstuff
baslineSV2(ampit)=mean(alldata{ampit}(sorted_indices(1:halflength),subselectbaseline)-mean(alldata{ampit}(sorted_indices(1:halflength),baselinetime),2,'omitnan'),'all','omitnan')*multiplyspk;
baselineSstdV2(ampit)=SEM(alldata{ampit}(sorted_indices(1:halflength),subselectbaseline)-mean(alldata{ampit}(sorted_indices(1:halflength),baselinetime),'all','omitnan'),0).*multiplyspk;
avgresp2V2(ampit)=mean(alldata{ampit}(sorted_indices(1:halflength),avgtime)-mean(alldata{ampit}(sorted_indices(1:halflength),baselinetime),2,'omitnan'),'all','omitnan')*multiplyspk;
stdresp2V2(ampit)=SEM(alldata{ampit}(sorted_indices(1:halflength),avgtime)-mean(alldata{ampit}(sorted_indices(1:halflength),baselinetime),2,'omitnan'),0).*multiplyspk;

baslineLV2(ampit)=mean(alldata{ampit}(sorted_indices(halflength:end),subselectbaseline)-mean(alldata{ampit}(sorted_indices(halflength:end),baselinetime),2,'omitnan'),'all','omitnan')*multiplyspk;
baselineLstdV2(ampit)=SEM(alldata{ampit}(sorted_indices(halflength:end),subselectbaseline)-mean(alldata{ampit}(sorted_indices(halflength:end),baselinetime),2,'omitnan'),0).*multiplyspk;
avgrespLV2(ampit)=mean(alldata{ampit}(sorted_indices(halflength:end),avgtime)-mean(alldata{ampit}(sorted_indices(halflength:end),baselinetime),2,'omitnan'),'all','omitnan')*multiplyspk;
stdrespLV2(ampit)=SEM(alldata{ampit}(sorted_indices(halflength:end),avgtime)-mean(alldata{ampit}(sorted_indices(halflength:end),baselinetime),2,'omitnan'),0).*multiplyspk;

%ttest(mean(alldata{ampit}(sorted_indices(halflength:end),avgtime)-mean(alldata{ampit}(sorted_indices(halflength:end),baselinetime),2,'omitnan')), mean(alldata{ampit}(sorted_indices(1:halflength),avgtime)-mean(alldata{ampit}(sorted_indices(1:halflength),baselinetime),2,'omitnan')))%test significance in V2 - V1 test on line 587

end

%sigmoidfig
avgtime=92:112;
figure(200)
hold on
baslineS=mean(HLB_data{ampit}{2}(:,subselectbaseline)-mean(HLB_data{ampit}{2}(:,baselinetime),2,'omitnan'),'all','omitnan')*multiplyspk;
baselineSstd=SEM(HLB_data{ampit}{2}(:,subselectbaseline)-mean(HLB_data{ampit}{2}(:,baselinetime),2,'omitnan'),0).*multiplyspk;
avgresp2=cellfun(@(x) mean(x{2}(:,avgtime)-mean(x{2}(:,baselinetime),2,'omitnan'),'all','omitnan')*multiplyspk,HLB_data);
stdresp2=cellfun(@(x) SEM(x{2}(:,avgtime)-mean(x{2}(:,baselinetime),2,'omitnan'),0)*multiplyspk,HLB_data);
errorbar([0; savefilename{15}{2}.AMP],[baslineS; avgresp2],[baselineSstd;stdresp2],'b')
baselineL=mean(HLB_data{ampit}{1}(:,subselectbaseline)-mean(HLB_data{ampit}{1}(:,baselinetime),2,'omitnan'),'all','omitnan').*multiplyspk;
baselineLstd=SEM(HLB_data{ampit}{1}(:,subselectbaseline)-mean(HLB_data{ampit}{1}(:,baselinetime),2,'omitnan'),0).*multiplyspk;
avgresp=cellfun(@(x) mean(x{1}(:,avgtime)-mean(x{1}(:,baselinetime),2,'omitnan'),'all','omitnan')*multiplyspk,HLB_data);
stdresp=cellfun(@(x) SEM(x{1}(:,avgtime)-mean(x{1}(:,baselinetime),2,'omitnan'),0)*multiplyspk,HLB_data);
errorbar([0; savefilename{15}{2}.AMP],[baselineL; avgresp],[baselineLstd; stdresp],'r')
set(gca,'TickDir','out');
xlabel('Current (\muA)');
ylabel('Average Firing Rate (Sp/s)')
title('V1')
legend('Bottom 1/2','Top 1/2')
%test signifiance - ttest(mean(HLB_data{ampit}{1}(:,avgtime)-baselineL,1,'omitnan'),mean((HLB_data{ampit}{2}(:,avgtime)-baslineS),1,'omitnan'))
figure(201)
hold on
errorbar([0; savefilename{15}{2}.AMP],[mean(baslineLV2); avgrespLV2],[mean(baselineLstdV2); stdrespLV2],'r')
errorbar([0; savefilename{15}{2}.AMP],[mean(baslineSV2); avgresp2V2],[mean(baselineSstdV2); stdresp2V2],'b')
set(gca,'TickDir','out');
xlabel('Current (\muA)');
ylabel('Average Firing Rate (Sp/s)')
title('V2')
legend('Bottom 1/2','Top 1/2')

%10ua timing figures
figure(100+ampit)
ax=axes;
hold on
baselineS=mean(alldata{ampit}(sorted_indices(1:halflength),1:89),2,'omitnan');
baselineL=mean(alldata{ampit}(sorted_indices(halflength+1:end),1:89),2,'omitnan');
stdshade((alldata{ampit}(sorted_indices(1:halflength),:)-baselineS).*multiplyspk,0.2,'b',[-90:90],1,ax);
xlabel('Time(ms)')
ylabel('FR (sp/s)')
title('V2')
 xlim([-50,89])
set(gca,'TickDir','out');
stdshade((alldata{ampit}(sorted_indices(halflength+1:end),:)-baselineL).*multiplyspk,0.2,'r',[-90:90],1,ax);
text(-40,3.5,['N=',num2str(halflength),'electrodes'],'Color','red')
text(-40,3.2,['N=',num2str(length(alldata{ampit}(sorted_indices(halflength+1:end),1))),'electrodes'],'Color','blue')
legend('Bottom 1/2','Top 1/2')
figure(105+ampit)
ax=axes;
hold on

baselineL=mean(HLB_data{ampit}{1}(:,1:89),2,'omitnan');
baselineS=mean(HLB_data{ampit}{2}(:,1:89),2,'omitnan');
stdshade((HLB_data{ampit}{1}-baselineL).*multiplyspk,0.2,'r',[-90:90],1,ax);
stdshade((HLB_data{ampit}{2}-baselineS).*multiplyspk,0.2,'b',[-90:90],1,ax);
%plot(-90:90,mean(HLB_data{3}.*multiplyspk,'omitnan'))
text(-40,20,['N=',num2str(size(HLB_data{ampit}{1},1)),'electrodes'],'Color','red')
text(-40,18,['N=',num2str(size(HLB_data{ampit}{2},1)),'electrodes'],'Color','blue')
xlabel('Time(ms)')
ylabel('FR (sp/s)')
title('V1')
 xlim([-50,89])
 ylim([-5 35])
set(gca,'TickDir','out');

%% early late
figure
ax=axes;
hold on
baselineE=mean(DataEarlyV2{5}(:,1:89),2);
baselineL=mean(DataLateV2{5}(:,1:89),2);
stdshade((DataEarlyV2{5}-baselineE),0.2,'r',[-90:90],1,ax);
stdshade((DataLateV2{5}-baselineL),0.2,'b',[-90:90],1,ax);
% plot(mean(DataEarlyV2{5}),'r')
% plot(mean(DataLateV2{5}),'b')
title('V2')
xlabel('Time (ms)')
ylabel('Firing rate (sp/s)')
 xlim([-50,89])
set(gca,'TickDir','out');
legend('Early','Late')


sigmoidELV1=nan(5,2);
sigmoidELstdV1=nan(5,2);
baselineV1EL=nan(5,2);
baselineELstdV1=nan(5,2);
sigmoidELV2=nan(5,2);
sigmoidELstdV2=nan(5,2);
baselineV2EL=nan(5,2);
baselineELstdV2=nan(5,2);
timesigmoid=92:122;
for ampit=1:5
    %V1
 dataE=reshape(mean(DataEarlyV1{ampit},2,'omitnan'),64,size(DataEarlyV1{ampit},1)/64);
 dataE(isnan(dataE))=0;
 dataL=reshape(mean(DataLateV1{ampit},2,'omitnan'),64,size(DataLateV1{ampit},1)/64);
 dataL(isnan(dataL))=0;
[~,iE] = setdiff(dataE.', dataL.', 'rows');%making sure its exclusively early or late ->NOT BOTH
[~,iL] = setdiff(dataL.', dataE.', 'rows');
dataE=[];
for i=1:length(iE)
    dataE=[dataE; DataEarlyV1{ampit}(iE(i)*64-63:iE(i)*64,:)];
end
dataE=unique(dataE,'rows');
dataL=[];
for i=1:length(iL)
    dataL=[dataL; DataLateV1{ampit}(iL(i)*64-63:iL(i)*64,:)];
end
dataL=unique(dataL,'rows');
baselineE=mean(dataE(:,baselinetime),2,'omitnan');
baselineL=mean(dataL(:,baselinetime),2,'omitnan');

sigmoidELV1(ampit,2)=mean(dataL(:,timesigmoid)-baselineL,'all','omitnan');
sigmoidELV1(ampit,1)=mean(dataE(:,timesigmoid)-baselineE,'all','omitnan');
sigmoidELstdV1(ampit,1)=SEM(dataE(:,timesigmoid)-baselineE,0);
sigmoidELstdV1(ampit,2)=SEM(dataL(:,timesigmoid)-baselineL,0);

baselineV1EL(ampit,2)=mean(dataL(:,subselectbaseline)-baselineL,'all','omitnan');
baselineV1EL(ampit,1)=mean(dataE(:,subselectbaseline)-baselineE,'all','omitnan');
baselineELstdV1(ampit,1)=SEM(dataE(:,subselectbaseline)-baselineE,0);
baselineELstdV1(ampit,2)=SEM(dataL(:,subselectbaseline)-baselineL,0);
%V2

baselineE=mean(DataEarlyV2{ampit}(:,baselinetime),2,'omitnan');
baselineL=mean(DataLateV2{ampit}(:,baselinetime),2,'omitnan');
sigmoidELV2(ampit,2)=mean(DataLateV2{ampit}(:,timesigmoid)-baselineL,'all','omitnan');
sigmoidELV2(ampit,1)=mean(DataEarlyV2{ampit}(:,timesigmoid)-baselineE,'all','omitnan');
sigmoidELstdV2(ampit,1)=SEM(DataEarlyV2{ampit}(:,timesigmoid)-baselineE,0);
sigmoidELstdV2(ampit,2)=SEM(DataLateV2{ampit}(:,timesigmoid)-baselineL,0);

baselineV2EL(ampit,2)=mean(DataLateV2{ampit}(:,subselectbaseline)-baselineL,'all','omitnan');
baselineV2EL(ampit,1)=mean(DataEarlyV2{ampit}(:,subselectbaseline)-baselineE,'all','omitnan');
baselineELstdV2(ampit,1)=SEM(DataEarlyV2{ampit}(:,subselectbaseline)-baselineE,0);
baselineELstdV2(ampit,2)=SEM(DataLateV2{ampit}(:,subselectbaseline)-baselineL,0);
end
figure
ax=axes;
hold on
baselineE=mean(dataE(:,1:89),2,'omitnan');
baselineL=mean(dataL(:,1:89),2,'omitnan');
stdshade((dataE-baselineE),0.2,'r',[-90:90],1,ax);
stdshade((dataL-baselineL),0.2,'b',[-90:90],1,ax);
% plot(mean(dataE,'omitnan'),'r')
% plot(mean(dataL,'omitnan'),'b')
title('V1')
xlabel('Time (ms)')
ylabel('Firing rate (sp/s)')
 xlim([-50,89])
set(gca,'TickDir','out');
legend('Early','Late')
ylim([-2.5 30])


%sigmoid early late
figure
hold on
errorbar([0 2 5 6 8 10],[mean(baselineV1EL(:,1)); sigmoidELV1(:,1)],[mean(baselineELstdV1(:,1)); sigmoidELstdV1(:,1)],'r')
errorbar([0 2 5 6 8 10],[mean(baselineV1EL(:,2)); sigmoidELV1(:,2)],[mean(baselineELstdV1(:,2)); sigmoidELstdV1(:,2)],'b')
set(gca,'TickDir','out');
legend('Early','Late')
title('V1')
xlabel('Current (ua)')
ylabel('Firing rate (sp/s)')
figure
hold on
errorbar([0 2 5 6 8 10],[mean(baselineV2EL(:,1)); sigmoidELV2(:,1)],[mean(baselineELstdV2(:,1)); sigmoidELstdV2(:,1)],'r')
errorbar([0 2 5 6 8 10],[mean(baselineV2EL(:,2)); sigmoidELV2(:,2)],[mean(baselineELstdV2(:,2)); sigmoidELstdV2(:,2)],'b')
set(gca,'TickDir','out');
legend('Early','Late')
title('V2')
xlabel('Current (ua)')
ylabel('Firing rate (sp/s)')
% baselineS{ampit}{1}=mean(HLB_data{ampit}{1}(:,1:89),2,'omitnan');
% baselineFR{ampit}{2}=mean(HLB_data{ampit}{2}(:,1:89),2,'omitnan');
% errorbar((savefilename{15}{2}.AMP(ampit)),mean(HLB_data{ampit}{1}(:,avgtime)-baselineFR{ampit}{1},'all','omitnan').*multiplyspk,'r')
% errorbar((savefilename{15}{2}.AMP(ampit)),mean(HLB_data{ampit}{2}(:,avgtime)-baselineFR{ampit}{2},'all','omitnan').*multiplyspk,'b')
%   
% correlation
% [maxbig,lagbig]=max(bigcorall{5}(:,181:end),[],2);
% lagbig(isnan(maxbig))=0;
% maxbig(isnan(maxbig))=0;
% lagbig((maxbig)==0)=[];
% maxbig((maxbig)==0)=[];
% lagbig((maxbig)==1)=[];
% maxbig((maxbig)==1)=[];
% figure
% scatter(lagbig,maxbig,'r')
% 
% [maxsmall,lagsmall]=max(smallcorall{5}(:,181:end),[],2);
% lagsmall(isnan(maxsmall))=0;
% maxsmall(isnan(maxsmall))=0;
% lagsmall((maxsmall)==0)=[];
% maxsmall((maxsmall)==0)=[];
% lagsmall((maxsmall)==1)=[];
% maxsmall((maxsmall)==1)=[];
% hold on
% scatter(lagsmall,maxsmall,'b')


% subsnum=1:length(lagsmall);
% permuted_numbers = subsnum(randperm(length(subsnum)));
% subselectnum=permuted_numbers(1:length(lagbig)); %baseline to subtract from

% 
% figure
% hold on
% 
% histogram(lagsmall(subselectnum))
% histogram(lagbig)
% title('lags')
% figure
% hold on
% histogram(maxsmall(subselectnum),8)
% histogram(maxbig,8)
% title('correlation')
% bigcorall
% smallcorall

%% early late
ampit=5;
figure
plot(mean(V2ELresp{ampit}{1}.*multiplyspk,'omitnan'),'r')
hold on
plot(mean(V2ELresp{ampit}{3}.*multiplyspk,'omitnan'),'b')

figure
plot(mean(V1ELresp{ampit}{1}.*multiplyspk,'omitnan'),'r')
hold on
plot(mean(V1ELresp{ampit}{3}.*multiplyspk,'omitnan'),'b')




%% sigmoid
baselintimee=1:89;%2:12 for V1, 1:89 for v2
peaktime=92:180;%92:102 for V1, 92:180 for v2
avgtimeV1=peaktime;
avgalldata=cellfun(@(x) mean(x(:,avgtimeV1).*multiplyspk,'all'),alldata);%avg of all time points
baselinedata=mean(cellfun(@(x) mean((x(:,baselintimee).*multiplyspk),'all'), alldata));
stdalldata=cellfun(@(x) SEM(x(:,avgtimeV1).*multiplyspk,0), alldata);
stdalldatabaseline=mean(cellfun(@(x) SEM(x(:,baselintimee).*multiplyspk,0), alldata));%
figure
%plot(savefilename{15}{2}.AMP,avgalldata)
errorbar([0;savefilename{15}{2}.AMP],[baselinedata;avgalldata],[stdalldatabaseline;stdalldata])
xlabel('Current (\muA)')
ylabel('Firing rate (sp/s)')
xlim([0 10])
set(gca,'TickDir','out');
hold on

avgalldata=cellfun(@(x) mean(max(x(:,peaktime).*multiplyspk,[],2)), alldata);%find peak of the data
baselinedata=mean(cellfun(@(x) mean(max(x(:,baselintimee).*multiplyspk,[],2)), alldata));
stdalldata=cellfun(@(x) SEM(max(x(:,peaktime).*multiplyspk,[],2),0), alldata);%
stdalldatabaseline=mean(cellfun(@(x) SEM(max(x(:,baselintimee).*multiplyspk,[],2),0), alldata));%
errorbar([0;savefilename{15}{2}.AMP],[baselinedata;avgalldata],[stdalldatabaseline;stdalldata])
legend('Avg','Peak')

%use this to compare ttest(mean(alldata{5}(:,avgtimeV1),1),mean(alldata{4}(:,avgtimeV1),1))

% avgalldata=cellfun(@(x) mean(x(:,peaktime).*multiplyspk,'all'),alldata);%avg of all time points
% %avgalldata=cellfun(@(x) mean(max(x(:,92:181).*multiplyspk,[],2)), alldata);%find peak of the data
% %baselinedata=mean(cellfun(@(x) mean(max(x(:,1:85).*multiplyspk,[],2)), alldata));
% baselinedata=mean(cellfun(@(x) mean((x(:,1:89).*multiplyspk),'all'), alldata));
% %stdalldata=cellfun(@(x) std(x(:,92:181).*multiplyspk,0,'all'), alldata);
% stdalldata=cellfun(@(x) std(max(x(:,peaktime).*multiplyspk,[],2),0,'all')./length(max(x(:,peaktime).*multiplyspk,[],2)), alldata);%
% stdalldatabaseline=mean(cellfun(@(x) std(max(x(:,baselintimee).*multiplyspk,[],2),0,'all')./length(max(x(:,baselintimee).*multiplyspk,[],2)), alldata));%


%% check where mean resp crosses each other in time
resp=cellfun(@(x) (mean(x.*multiplyspk)),alldata,'UniformOutput',false);
respall=[resp{1}-resp{2};resp{1}-resp{3};resp{1}-resp{4};resp{1}-resp{5}];
[~,c]=find(respall(4,94:181)>0,5,'first');

%% stimchnalilgn heatmap
datastimchn=mean(centredstimchn{5}.*1000,3,'omitnan');
figure
imagesc(datastimchn)
colorbar

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

%% Distance heatmap response
figure; heatmap(sum(heatmap_centroid{5}>0,3)); set(gca,'ColorScaling','log')
dataheatmap=sum(heatmap_centroid{5}>0,3);
figure; heatmap(dataheatmap,'CellLabelColor','none','GridVisible','off');set(gca,'ColorScaling','log')
maxdata=max(dataheatmap,[],'all')/2;
fwhmdat=dataheatmap>=maxdata;

%stimchn only
dataheatmap=sum(centredstimchn{5}>0,3);
figure; heatmap(dataheatmap,'CellLabelColor','none','GridVisible','off');set(gca,'ColorScaling','log')


numRows = 31;  % Number of rows (each row is 50 units apart)
numCols = 7;   % Number of columns (each column is 200 units)
% Create an array of row indices
rowIndices = (-floor(numRows/2):floor(numRows/2)) * 50;
% Create an array of column indices
colIndices = (-floor(numCols/2):floor(numCols/2)) * 200;
% Create the distance array
distanceArray = sqrt((rowIndices.').^2 + (colIndices).^2);
maxnum=1000;
step=100;
splitdist=cell(5,1);
splitdist(cellfun(@isempty, splitdist)) =   {nan((maxnum/step),1)};
stdsplitdist=splitdist;
countsplidist=splitdist;
figure; hold on;
if all(chnrange>64)
dataplot=centredstimchn;
xlabel('Distance from stim electrode(\mum)')
elseif all(chnrange<65)
    dataplot=heatmap_centroid;
    xlabel('Distance from centroid(\mum)')
end
    
for current=1:2:5
    dataheatmap=sum(dataplot{current}>0,3,'omitnan');
for i=step:step:maxnum-step
    splitdist{current}(i/step)=mean(dataheatmap(distanceArray>=i-step & distanceArray<i));
    stdsplitdist{current}(i/step)=std(dataheatmap(distanceArray>=i-step & distanceArray<i));
    countsplidist{current}(i/step)=sum(dataheatmap(distanceArray>=i-step & distanceArray<i));
   
end
 errorbar(step:step:maxnum,splitdist{current},stdsplitdist{current})
end
color1 = linspace(0,1,3);
newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
colororder(newcolors);
set(gca,'TickDir','out');

ylabel('# significant electrodes')
% %% sigmoid using groupdata
% sigmoidampepoch=zeros(5,9);
% for ampit=1:5
%     for epoch=1:9
%         sigmoidampepoch(ampit,epoch)=mean(alldata10uA{ampit}{epoch}(:,82+(10*epoch):91+(10*epoch)),'all'); 
%     end
% end
% figure; hold on
% for i=1:9
% plot([2 5 6 8 10],sigmoidampepoch(:,i).*multiplyspk)
% end
% color1 = linspace(0,1,9);
% newcolors = [zeros(length(color1),1) flipud(color1') (color1')];
% colororder(newcolors);
% %plot([2 5 6 8 10],mean(sigmoidampepoch,2).*1000,'k')
% xlabel('current (\muA)')
% ylabel('Firing rate (sp/s)')
% leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','61:72','72:81','82:91','Average');
% %leg=legend('2:11','22:31','42:51','61:72','82:91');
% title(leg,'Time epoch(ms)')
% set(gca,'TickDir','out');
% 
% indvsig=zeros(5,9);
% for i=1:5
% indvsig(i,:)=cellfun(@(x) size(x,1), groupdata{i});
% end
% 
% figure
% hold on
% plot([2 5 6 8 10],indvsig)
% colororder(newcolors);
% plot([2 5 6 8 10],sum(indvsig,2),'k')
% yline((0.003)*totalchncount(1),'r')
% xlabel('current (\muA)')
% ylabel('# elect sig')
% %leg=legend('2:11','22:31','42:51','62:71','82:91','Total','Significance');
% leg=legend('2:11','12:21','22:31','32:41','42:51','52:61','61:72','72:81','82:91','Total','Significance');
% title(leg,'Time epoch(ms)')
% set(gca,'TickDir','out');

end
