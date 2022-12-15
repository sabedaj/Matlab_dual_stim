%%brain stimulation figure generation - original code was
%%loopsubdirectpoolactivity

%% for rasterplots
%flash
%E:\DATA\Rat_019\flash_001_210511_100735 chn 9 original
%E:\DATA\Rat_014\S1E2_9elect_001_210311_165352 chn 29
count_sp_tr=Flash_raster(29) %only plotted the first 200 trials, this will plot all 900


%stim
%E:\DATA\Rat_019\S1E4_9elect_h_001_210511_161244 chn 9
%E:\DATA\Rat_014\S1E2_9elect_001_210311_165352 chn 29
%generate_StackedRaster_sab(9,6) %plotted chn 9, 6uA original
generate_StackedRaster_sab(29,6) %plotted chn 9, 6uA
% 
% for i=26:32
% generate_StackedRaster_sab(i,6)
% end
% plotting individual spike traces
% middle three E:\DATA\Rat_019\S1E4_9elect_h_001_210511_161244 tID==13
% top %E:\DATA\Rat_019\S1E5_7elect_h_001_210511_184045 tID>5
% bottom E:\DATA\Rat_019\S3E2_9elect_001_210511_104603 tID>5
Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
startpointseconds=2; %How long after the trigger do you want skip spike analysis(ms)? 
secondstoanalyse=8; %How long after the trigger do you want to analyse spikes for(ms)? 
printspiking=0;
trig = loadTrig(0);
TrialParams=loadTrialParams;
maxid=max(cell2mat(TrialParams(:,2)));
starttrial=1;
trialjump=1;
endtrial=maxid;
%for spike traces
%tID==13 && find(E_MAP==chsp)==50   E:\DATA\Rat_019\S1E4_9elect_h_001_210511_161244
[IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);

%% Initial loop through data to organise it for plotting
%for half amplitude, change poolnormalised activity flag at beginning of
%function. To change to baseline activity, alter Csplit line 216 to basespk_all
%instead of normspk_all. Single animal, just change ratN to the animal - R16
%for 5 elect and R19 for 7,9 elect.
marmoflag=0;
startpointseconds=2;
secondstoanalyse=8;
numdistsigelect_array=NaN(21,1700);
distsigelect_array=NaN(21,1700);
spkrate_dist=NaN(21,1700);
counterstimchn=0;
counternumberofpairs=0;
stimpos=[];
for ratN=14:20
    Ratnum=['Rat_0' num2str(ratN)];
saveshank.(Ratnum)=[];
stimchn_r.(Ratnum)=[];
end
countchn=0;
countchn(2)=0;
countloop=0;
for sepdist=-6:-2:-10 %%% Laminar
AMP=[0 1 2 3 4 5 6 8 10 15];
sepcheck=['sep' num2str((sepdist+1)*-1)];
rollingsum.(sepcheck)=nan(16,1000);
altogethersumforavg.(sepcheck)=nan(1,1000);
for trial=1:5
    check1=['T' num2str(trial)];
    for current = 1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        Csplit.(currcheck).(check1)=[];
        basesplit.(currcheck).(check1)=[];
        sepdistsig_dual.(sepcheck).(currcheck)=nan(32,1000);
        sepdistsig.(sepcheck).(currcheck)=nan(32,1000);
        CsplitCentroid.(currcheck)=[];
        for shanksep=0:3
            bincheck=['D' num2str(shanksep)];
            bin.(currcheck).(check1).(bincheck)=[];
%             for TimeBegin=2:2:6
%                 timecheck=['T' num2str(TimeBegin)];
%                 bin_time.(currcheck).(check1).(bincheck).(timecheck)=[];
%             end
        end
    end
end

countershank=0;
stimshank=[];
AR=[1 2 3 4];%array shanks
binchecksave='blank';

for ratN=14:20
if ratN<10 && marmoflag==0
    Ratnum=['Rat_00' num2str(ratN)];
elseif ratN>=10 && marmoflag==0
    Ratnum=['Rat_0' num2str(ratN)];
else%%marmo
    Ratnum=['Marmo_00' num2str(ratN)];
end
AMP=[0 1 2 3 4 5 6 8 10 15];

cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
D_data=dir;

savestimshanks=zeros(length(D_data),1);
       
for k = 3:length(D_data) % avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    try
        cd([D_data(k).folder filesep currD])
        loadStimChn;
        if (stimChn(1)-stimChn(2))~=sepdist
            continue
        end
        load('Averagetrialresponse.mat','avgnospT')
    catch
        stop=0;
        continue
    end

        [normspk_all,centroidpershank, basespk_all, respspk_all]=PoolNormalisedActivity(avgnospT,startpointseconds, secondstoanalyse);
        if sum(normspk_all.S1_T1,'all')==0 && sum(centroidpershank,'all')==140
            continue
        end

        stimChn_orig=stimChn;
        if stimChn(1)<17 %determines the shank with the stimulating electrodes
            shank=1;
        elseif stimChn(1)<33 && stimChn(1)>16
            shank=4;
            stimChn=stimChn-16;
        elseif stimChn(1)<49 && stimChn(1)>32
            shank=2;
            stimChn=stimChn-32;
        else
            shank=3;
            stimChn=stimChn-48;
        end
        [~,uniquestimchn]=intersect((stimChn_orig),stimchn_r.(Ratnum));%multiply by shanknum to ensure unique single stimulation channel
        countloop=countloop+1;
        if (isempty(uniquestimchn))
            countchn(1)=countchn(1)+2;
        elseif length(uniquestimchn)==1
            countchn(1)=countchn(1)+1;
            countchn(2)=countchn(2)+1;
        else
            countchn(2)=countchn(2)+2;
        end
        %stimpos=[stimpos; stimChn];
        AMP=loadAMP;
        AMP(AMP==-1)=0;

            clear baslinespikestruct IDstruct
            load('IDstruct.mat')
            Ampcheck=loadAMP;
            if any(Ampcheck==6) && exist('baslinespikestruct')==1
                counterstimchn=counterstimchn+1;
            end
            trialinfo=loadTrialInfo(0);
            trialinfo(:,3)=[];
            trialinfo=cell2mat(trialinfo);
         for shanknum=1:4
                if shanknum==2 %used to index data as it is stored shank 1,4,2,3
                    shank_orig=3;
                elseif shanknum==3
                    shank_orig=4;
                elseif shanknum==4
                    shank_orig=2;
                else
                    shank_orig=1;
                end
            countershank=countershank+1;
            trialcheck=0;
                        
            for trial=1:5
                check=['S' num2str(shanknum) '_T' num2str(trial)];
                % Identifying the number of significantly responding electrodes based on shank, distance from stimulating electrode, dual and single electrode conditions
                if (trial==1 ||trial==5)%&& (isempty(uniquestimchn) || all(uniquestimchn~=1)))  || (trial==5 && (isempty(uniquestimchn) || all(uniquestimchn~=2))) %single electrode trials only
                    if (trial==5)
                        chnstim=2;
                    else
                        chnstim=1;
                    end
                    if any(Ampcheck==6) && exist('baslinespikestruct')==1
                        trialsAmpmatch=find((trialinfo(:,2)==stimChn_orig(chnstim)) & (trialinfo(:,18)==6));
                        [trialsAmpmatchindex]=find(trialinfo(trialsAmpmatch+1,2)==0);
                        trialmatch=trialinfo(trialsAmpmatch(trialsAmpmatchindex),1);
                        ID=trialmatch;
                        checkID=['T' num2str(ID)];
                        E_MAP = Depth(1);
                        IDstructID=IDstruct.(checkID)(E_MAP,:);
                        baslinespikestructID=baslinespikestruct.(checkID)(E_MAP,:); %comparing baseline to stimualtion trials to determine significance
                        for numelect=1:size(basespk_all.(check),1)
                            [psigchan,hsigchan,~] = signrank(baslinespikestructID(numelect*shank_orig,:), IDstructID(numelect*shank_orig,:),'alpha',0.05,'tail','left');
                            rollingsum.(sepcheck)(numelect,(counterstimchn*2)-(chnstim-1))=nansum([rollingsum.(sepcheck)(numelect,(counterstimchn*2)-(chnstim-1)) mean(IDstructID(numelect*shank_orig,:))*1000/6]);%% all significant channels %hsigchan]); %
                            distancerecordstim=round(sqrt(((stimChn(chnstim)-numelect)*50)^2+((shank-shanknum)*200)^2)/50); %distance of the electrode from the single stimulating electrode sorted into 50um bins
                            numdistsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1))=nansum([numdistsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1)) 1]);
                            distsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1))=nansum([distsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1)) mean(IDstructID(numelect*shank_orig,:))*1000/6]);% hsigchan]);%
                            if hsigchan==1
                                spkrate_dist(distancerecordstim+1,(counterstimchn*2)-(chnstim-1))=nansum([spkrate_dist(distancerecordstim+1,(counterstimchn*2)-(chnstim-1)) normspk_all.(check)(numelect,Ampcheck==6)]);
                            end
                        end

                        altogethersumforavg.(sepcheck)((counterstimchn*2)-(chnstim-1))=nansum([altogethersumforavg.(sepcheck)((counterstimchn*2)-(chnstim-1)) 1]);
                    elseif any(Ampcheck==6)
                        stop=0;
                    end
                    if exist('baslinespikestruct')==1
                        for current = 1:length(AMP)
                            if AMP(current)==0 || ~any(Ampcheck==AMP(current))
                                continue
                            end
                            currcheck=['C' num2str(AMP(current))];
                            trialsAmpmatch=find((trialinfo(:,2)==stimChn_orig(chnstim)) & (trialinfo(:,18)==AMP(current)));
                            [trialsAmpmatchindex]=find(trialinfo(trialsAmpmatch+1,2)==0);
                            trialmatch=trialinfo(trialsAmpmatch(trialsAmpmatchindex),1);
                            ID=trialmatch;
                            checkID=['T' num2str(ID)];
                            E_MAP = Depth(1);
                            IDstructID=IDstruct.(checkID)(E_MAP,:);
                            baslinespikestructID=baslinespikestruct.(checkID)(E_MAP,:);
                            for numelect=1:size(basespk_all.(check),1)
                                [psigchan,hsigchan,~] = signrank(baslinespikestructID(numelect*shank_orig,:), IDstructID(numelect*shank_orig,:),'alpha',0.05,'tail','left'); 
                                sepdistsig.(sepcheck).(currcheck)(numelect+(16-stimChn(1)),(counterstimchn*2)-(chnstim-1))=nansum([sepdistsig.(sepcheck).(currcheck)(numelect+(16-stimChn(1)),(counterstimchn*2)-(chnstim-1)) hsigchan]);
                            end
                        end
                    end

                    
                elseif ((isempty(uniquestimchn) || all(uniquestimchn~=1)) || (isempty(uniquestimchn)|| all(uniquestimchn~=2))) && (trial~=1 && trial~=5)
                    if exist('baslinespikestruct')==1
                         for current = 1:length(AMP)
                             if AMP(current)==0 || ~any(Ampcheck==AMP(current))
                                 continue
                             end
                             currcheck=['C' num2str(AMP(current))];
                            trialsAmpmatch=find((trialinfo(:,2)==stimChn_orig(1)) & (trialinfo(:,18)==(trial-1)*0.25*AMP(current)))+1;
                            trialsAmpmatch2=find((trialinfo(:,2)==stimChn_orig(2)) & (trialinfo(:,18)==(5-trial)*0.25*AMP(current)));
                            trialsAmpmatchindex=intersect(trialsAmpmatch,trialsAmpmatch2);
                            trialmatch=trialinfo(trialsAmpmatchindex,1);
                            ID=trialmatch;
                            checkID=['T' num2str(ID)];
                            E_MAP = Depth(1);
                            IDstructID=IDstruct.(checkID)(E_MAP,:);
                            baslinespikestructID=baslinespikestruct.(checkID)(E_MAP,:);
                            for numelect=1:size(basespk_all.(check),1)
                                [psigchan,hsigchan,~] = signrank(baslinespikestructID(numelect*shank_orig,:), IDstructID(numelect*shank_orig,:),'alpha',0.05,'tail','left');
                                chnstim=1;
                                sepdistsig_dual.(sepcheck).(currcheck)(numelect+(16-stimChn(chnstim)),(counterstimchn*3)-(trial-2))=nansum([sepdistsig_dual.(sepcheck).(currcheck)(numelect+(16-stimChn(chnstim)),(counterstimchn*3)-(trial-2)) hsigchan]);
                            end
                         end
                    end
                end
                
                %divide up into time 2ms
%                 TrialParams=loadTrialParams;
%                 endtrial=max(cell2mat(TrialParams(:,2)));
%                 trig = loadTrig(0);
%                 TimeStep=2;
%                 for TimeBegin=2:2:6
%                 [IDstruct, baslinespikestruct]=sortTrials_SM(TimeBegin,TimeStep+TimeBegin,trig,0,1,1,endtrial,0);
%                 [avgnospT,~,~] = AverageTrialResponse_SM(IDstruct, baslinespikestruct);
%                 [normspk_all_time,~, ~, ~]=PoolNormalisedActivity(avgnospT,TimeBegin, TimeStep+TimeBegin);
%                 timecheck=['T' num2str(TimeBegin)];
%                 normspk_all_timesort.(timecheck)=normspk_all_time;
%                 end
                % spike rate and centroid calculation
                check1=['T' num2str(trial)];
                for current = 1:length(AMP)
                    currcheck=['C' num2str(AMP(current))];
                    basesplit.(currcheck).(check1)=[basesplit.(currcheck).(check1) [NaN(16-stimChn(1),1); basespk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                    Csplit.(currcheck).(check1)=[Csplit.(currcheck).(check1) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                    if trialcheck~=1
                        CsplitCentroid.(currcheck)=[CsplitCentroid.(currcheck) centroidpershank((shanknum-1)*5+1:(shanknum*5),current)-stimChn(1)];
                    end
                     % plot resutls of 200um away, 400um away and 600um away
                        Shanksep=shanknum-shank;
                        bincheck=['D' num2str(abs(Shanksep))];
                        bin.(currcheck).(check1).(bincheck)=[bin.(currcheck).(check1).(bincheck) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
%                         for TimeBegin=2:2:6
%                             timecheck=['T' num2str(TimeBegin)];
%                             bin_time.(currcheck).(check1).(bincheck).(timecheck)=[bin_time.(currcheck).(check1).(bincheck).(timecheck) [NaN(16-stimChn(1),1); normspk_all_time.(timecheck).(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
%                         end
%                         if strcmp(bincheck,binchecksave)
%                             bin.(currcheck).(check1).(bincheck)(:,end)=mean([bin.(currcheck).(check1).(bincheck)(:,end) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]],2);
%                         else
%                             bin.(currcheck).(check1).(bincheck)=[bin.(currcheck).(check1).(bincheck) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
%                         end
                end
                trialcheck=1;
            end
            if ~strcmp(bincheck,'D0')
                binchecksave=bincheck;
            end
         end
         binchecksave='blank';
         counternumberofpairs=counternumberofpairs+1;
         stimpos=[stimpos; stimChn];
        countershank=0;
        
        stimchn_r.(Ratnum)=[stimchn_r.(Ratnum) (stimChn_orig)];
        
        saveshank.(Ratnum)=[saveshank.(Ratnum) shank];
end
savestimshanks(savestimshanks==0)=[];
stimshank=[stimshank; savestimshanks];
end
check=['sep' num2str(-1*sepdist-1)];
saveCplit.(check)=Csplit;
savebasesplit.(check)=basesplit;
saveCsplitcentroid.(check)=CsplitCentroid;
saveshanksepdist.(check)=bin;
end
%% number of times a shank was stimulated - all stim electrodes not just unique
for ratN=14:20
    Ratnum=['Rat_0' num2str(ratN)];
 c = unique(saveshank.(Ratnum)); % the unique values in the A (1,2,3,4,5)    
 for i = 1:length(c)
   counts(i,ratN-13) = sum(saveshank.(Ratnum)==c(i)); % number of times each unique value is repeated
 end
end
mean(counts,'all')
std(counts,0,'all')

%pairs per animal
mean(sum(counts,1))
std(sum(counts,1))
%% histogram stim elect - all stim electrodes, not just unique
histogram(stimpos)
set(gca,'TickDir','out');
xlim([0 16])
ylabel('# stimulating electrodes')
xlabel('Distance from tip electrode (/mum)')
%% calculating the number of significantly responding electrodes
%57 stim pairs @ 6uA, 114 electrodes, not all unique!
allresp=[rollingsum.sep5 rollingsum.sep7 rollingsum.sep9];
alltog=[altogethersumforavg.sep5 altogethersumforavg.sep7 altogethersumforavg.sep9];
figure(3)
hold on
stdshade((allresp'./alltog'),0.2,'k',0:50:15*50)
%ylim([0 100])
xlim([0 750])
text(10,90,['N=' num2str(nansum(alltog)) ' shanks'])
%ylabel('Precentage responding electrodes (%)')
ylabel('Firing rate (sp/s)')
xlabel('Distance from tip electrode (um)')
title('Electrodes responding to 6uA stim compared to baseline')
set(gca,'TickDir','out');

% average num responding electrodes
nanmean((allresp./alltog).*100,'all')
nanstd(sum((allresp),1))/sqrt(sum(~isnan(sum((allresp),1))))

distancefrelect=0:50:18*50;
figure(4)%distance based sig elect
hold on
stdshade([distsigelect_array(1:19,:)./numdistsigelect_array(1:19,:)]',0.2,'k',distancefrelect)
%ylim([0 100])
%ylabel('Precentage responding electrodes (%)')
ylabel('Firing rate (sp/s)')
xlabel('Distance from stimulating electrode (\mum)')
set(gca,'TickDir','out');
title('Electrodes responding to 6uA stim compared to baseline')

%% Plotting fig 3 - main result -all shanks

singleCurrent=6; %plot 6uA results
ploty=1; %Do the plotting
if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end
vec = [100;80;50;30;15;0];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
centroidpos_all=[];
samples_rand=[];
singleelectspread=[];
dualelectspread=[];
trials=[1 2 3 4 5];
clear spread_electno spread_electnomaxall sepdistsigdwnsampe sepdistsigdwnsampe_dual peak_all
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    dualelectspreadsepdist.(sepcheck)=[];
    Csplit=saveCplit.(sepcheck);
    hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
    raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
    map = interp1(vec,raw,linspace(100,0,N),'pchip');
    cmap=colormap(map);
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        if length(AMP)>1 && ploty==1
            figure(1+sepdist)
            subplot(2,4,current)
            hold on
        end
        savdwnsamplespkrate=zeros(32,5);
        if ploty==1 && length(AMP)==1
            figure(1+sepdist)
            axes('Position',[0.13         0.112396822842341                     0.775         0.62])
            hold on
            set(gca,'FontSize',12)
        end
        for trial=trials
            check=['T' num2str(trial)];
            dat=Csplit.(currcheck).(check);
            dat(isinf(dat))=nan;
            clear pairavgcur erpairavgcur
            for paircount=1:size(dat,2)/4
                pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
            end
            lengthneeded=size(Csplit.C1.T1,2)/4;% find the smallest number of pairs for any current - we need to down-sample the others to match
            if trial==1 
                [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
                samples_rand.(sepcheck)=samples_rand1;
            end
            pairavg_dwnsample=pairavgcur(:, samples_rand.(sepcheck));% downsample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan; %remove rows with less than 10 pairs in the average
            pairavgdwnsample_all_nomean.(sepcheck).(currcheck).(check)=pairavg_dwnsample;
            pairavgdwnsample_all.(sepcheck).(currcheck)(:,trial)=nanmean(pairavg_dwnsample,2);
            if ploty==1
                stdshade(pairavg_dwnsample',0.2,cmap(trial*floor((length(cmap))/5),:));
            end
            check=['T' num2str(trial)];
         end

        
        if ploty==1
            ylim([0, round(max(pairavgdwnsample_all.(sepcheck).(currcheck),[],'all')/50)*50+50])
            if sepdist==5
                 xlim([16-3 16+sepdist+1+3])
            elseif sepdist==7
                 xlim([16-2 16+sepdist+1+2])
            elseif sepdist==9
                 xlim([16-1 16+sepdist+1+1])
            end
            xline(16,'r')
            xline(16+sepdist+1,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            yt = yticks;
            yticklabels(yt(1:end-1));
            ylabel('Sp/s')
            xlabel('Distance from deepest stim elect (\mum)')
            set(gca,'TickDir','out');
            ax1=gca;
            legend('0:100','25:75','50:50','75:25','100:0')
        end
        
    end
    if ploty==1
        set(gca,'TickDir','out');
        hex = ['#1b0c36';'#532e9e';'#147df5';'#02a612';'#fffb7d';'#ffe26e'];
        raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
        map = interp1(vec,raw,linspace(100,0,N),'pchip');
        cmap=colormap(map);
    end
    
    
    if (exist('sepdistsig'))==1 %checks spread and significance of single vs dual
        dualspread=[];
        singlespread=[];
        figure (30)
        subplot(1,3,round((sepdist-4)/2))
        hold on
        for paircount=1:size(sepdistsig.(sepcheck).C6,2)/2
            sepdistsigdwnsampe.(sepcheck)(:,paircount)=nanmean(nanmean(sepdistsig.(sepcheck).C6(:,( paircount-1)*2+1:2* paircount),2));
        end
        for paircount=1:size(sepdistsig_dual.(sepcheck).C6,2)/3
           sepdistsigdwnsampe_dual.(sepcheck)(:,paircount)=nanmean(nanmean(sepdistsig_dual.(sepcheck).C6(:,( paircount-1)*3+1:3* paircount),2));
        end
        sepdistsigdwnsampe_dual.(sepcheck)=[sepdistsigdwnsampe_dual.(sepcheck) nan(1, size(sepdistsigdwnsampe.(sepcheck),2)-size(sepdistsigdwnsampe_dual.(sepcheck),2))];
        psepsig.(sepcheck)=signrank(sepdistsigdwnsampe_dual.(sepcheck),sepdistsigdwnsampe.(sepcheck),'tail','right');
        ratiovar=sepdistsigdwnsampe_dual.(sepcheck)./sepdistsigdwnsampe.(sepcheck);
        ratiovar(isinf(ratiovar))=nan;
        avgratio.(sepcheck)=nanmean(ratiovar,'all')*100-100;
    end
    
    % centroid per pair %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcur=[];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];

        for trial=trials
            check=['T' num2str(trial)];
            dat=Csplit.(currcheck).(check);
            dat(isinf(dat))=nan;
            dat(sum(~isnan(dat),2)<3,:)=nan;
            clear pairavgcur erpairavgcur
            for paircount=1:size(dat,2)/4
                pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*4+1:4* paircount),2);
                [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur(:,paircount),sepdist);
                centroidpos_all.(sepcheck).(currcheck)(trial,paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%
            end
            lengthneeded=size(Csplit.C1.T1,2)/4;
            if trial==1
                [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
                samples_rand.(sepcheck)=samples_rand1;
            end
        end
        centroidpos_all.(sepcheck).(currcheck)=centroidpos_all.(sepcheck).(currcheck)(:, samples_rand.(sepcheck));
        centroidpos_all.(sepcheck).(currcheck)(sum(~isnan(centroidpos_all.(sepcheck).(currcheck)),2)<10,:)=nan;
        avgpcurr=mean(centroidpos_all.(sepcheck).(currcheck),2);
        stdavg=std(centroidpos_all.(sepcheck).(currcheck),[],2)./sqrt(size(centroidpos_all.(sepcheck).(currcheck),2));

        if ploty==1
            avgpcurr=avgpcurr(trials);
            if length(AMP)==1
            figure(1+sepdist)
            axes('Position',[ax1.Position(1) .73 ax1.Position(3) .2])
            ax2=gca;
            hold on
            set(gca,'FontSize',12)
            for colorbar=1:length(avgpcurr)
                color_current=cmap(trials(colorbar)*floor((length(cmap))/5),:);
                er = errorbar(avgpcurr(colorbar),colorbar,stdavg(trials(colorbar)),stdavg(trials(colorbar)),'horizontal');er.Color = [0, 0, 0]; er.LineStyle = 'none';
                scatter(avgpcurr(colorbar), colorbar, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)
            end
            yticks([0 1 2 3 4 5 6 7 8])
            labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
            labels=labels(trials);
            yticklabels([{''} labels {''}])
            ylim([0.5 length(trials)+0.5])
            xlim((ax1.XLim-16).*50)
            set(gca,'TickDir','out');
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            xline(0,'r')
            xline((sepdist+1)*50,'r')
            end
        end
    end
end


%% current vs sepdist
AMP=[0 1 2 3 4 6 8 10];
colorchoice={'r' 'b' 'k'};
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=2:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        currsepavgshift.(sepcheck) (current)=abs(mean(diff(peak_all.(sepcheck).(currcheck)),'all'));
        currseperrorshift.(sepcheck) (current)=std(diff(peak_all.(sepcheck).(currcheck)),[],'all')./sqrt(numel(diff(peak_all.(sepcheck).(currcheck))));
    end
    figure(5)
    hold on
    errorbar(AMP(2:end),currsepavgshift.(sepcheck)(2:end),currseperrorshift.(sepcheck)(2:end), colorchoice{sepdist/2-1.5})
end
ylabel('Avg centroid shift (um)')
xlabel('Current (uA)')
legend('300','400','500')
xlim([0 11])
ylim([0 100])
set(gca,'TickDir','out');
%% current vs spread
clear sepdistsigdwnsampe_dual sepdistsigdwnsampe
seedpoint=65;%random but needs to be consistent
s = RandStream('mlfg6331_64','Seed',seedpoint);
AMP=[0 1 2 3 4 6 8 10];
colorchoice={'r' 'b' 'k'};
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];

    for current=2:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for paircount=1:size(sepdistsig_dual.(sepcheck).(currcheck),2)/3
           sepdistsigdwnsampe_dual.(sepcheck).(currcheck)(:,paircount)=nanmean(nanmean(sepdistsig_dual.(sepcheck).(currcheck)(:,( paircount-1)*3+1:3* paircount),2));
        end
        sepdistsigdwnsampe_dual.(sepcheck).(currcheck)(isnan(sepdistsigdwnsampe_dual.(sepcheck).(currcheck)))=[];
        lengthneeded=length(sepdistsigdwnsampe_dual.(sepcheck).C1);
        [samples_rand1]=DownSample(sepdistsigdwnsampe_dual.(sepcheck).(currcheck),lengthneeded,s,seedpoint);
        samples_rand.(sepcheck)=samples_rand1;
        sepdistsigdwnsampe_dual.(sepcheck).(currcheck)=sepdistsigdwnsampe_dual.(sepcheck).(currcheck)(samples_rand1);
        totalpercentage_dualspread.(sepcheck)(current-1)=mean(sepdistsigdwnsampe_dual.(sepcheck).(currcheck))./4.*100;
        totalpercentagerrore_dualspread.(sepcheck)(current-1)=std(sepdistsigdwnsampe_dual.(sepcheck).(currcheck)./4.*100)./sqrt(numel(sepdistsigdwnsampe_dual.(sepcheck).(currcheck)));
        for paircount=1:size(sepdistsig.(sepcheck).(currcheck),2)/2
            sepdistsigdwnsampe.(sepcheck).(currcheck)(:,paircount)=nanmean(nanmean(sepdistsig.(sepcheck).(currcheck)(:,( paircount-1)*2+1:2* paircount),2));
        end
        sepdistsigdwnsampe.(sepcheck).(currcheck)(isnan(sepdistsigdwnsampe.(sepcheck).(currcheck)))=[];
        sepdistsigdwnsampe.(sepcheck).(currcheck)=sepdistsigdwnsampe.(sepcheck).(currcheck)(samples_rand1);
        totalpercentage_singlespread.(sepcheck)(current-1)=mean(sepdistsigdwnsampe.(sepcheck).(currcheck))./4.*100;
         totalerrorpercentage_singlespread.(sepcheck)(current-1)=std(sepdistsigdwnsampe.(sepcheck).(currcheck)./4.*100)./sqrt(numel(sepdistsigdwnsampe.(sepcheck).(currcheck)));
    end
    figure(6)
    hold on
     errorbar(AMP(2:end),totalpercentage_dualspread.(sepcheck),totalpercentagerrore_dualspread.(sepcheck), colorchoice{sepdist/2-1.5})
    h= errorbar(AMP(2:end),totalpercentage_singlespread.(sepcheck),totalerrorpercentage_singlespread.(sepcheck), ['--' colorchoice{sepdist/2-1.5}]);
    hold off
    alph=0.5;
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alph])
end
ylabel('Percentage electrodes responding to stim')
xlabel('Current')
legend('300 dual','300 single','400 dual','400 single','500 dual','500 single')
ylim([0 100])
xlim([0 11])
set(gca,'TickDir','out');
%% shift as a percentage of sep dist

%centroid
% datasep5=centroidpos_all.sep5.C6;
% datasep7=centroidpos_all.sep7.C6;
% datasep9=centroidpos_all.sep9.C6;

% %peak
datasep5=peak_all.sep5.C6;
datasep7=peak_all.sep7.C6;
datasep9=peak_all.sep9.C6;
trialstoplot=1:5;
figure
hold on
d1=(flipud(datasep5(trialstoplot,:)-150)'./(150).*100+100)./2;
d2=(flipud(datasep7(trialstoplot,:)-200)'./(200).*100+100)./2;
d3=(flipud(datasep9(trialstoplot,:)-250)'./(250).*100+100)./2;
stdshade(d1,0.2,'r',trialstoplot);
stdshade(d2,0.2,'b',trialstoplot);
stdshade(d3,0.2,'k',trialstoplot);

x=trialstoplot;
y=mean([d1; d2; d3],1);
p = polyfit(x,y,1); 
yfit = polyval(p,x); 
plot(x,yfit,'g')
numparam=2;
[rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,y,numparam);

xticks(1:5)
labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
labels=labels(trialstoplot);
xticklabels(labels)
ylabel('Peak location')
xlabel('% current delivered')
ylim([0 100])
set(gca,'TickDir','out');
legend('300\mum', '400\mum', '500\mum', 'Fit')
labels2={'Stim elect (tip)', 'Midway', 'Stim elect (base)'};
yticks([0 50 100])
yticklabels(labels2)


figure
hold on
stdshade(flipud(datasep5-150)',0.2,'r',1:5);
stdshade(flipud(datasep7-200)',0.2,'b',1:5);
stdshade(flipud(datasep9-250)',0.2,'k',1:5);
xticks(1:5)
labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
labels=labels(trials);
xticklabels(labels)
ylabel('Centroid position (\mum)')
xlabel('% current delivered')
set(gca,'TickDir','out');
legend('300\mum', '400\mum', '500\mum')

%% Plotting fig 3 - main result - single shank

plotAllShanks=1; %1=yes, 0=no. 0 only plots shanks directly next to stim shank
numshanksToAvg=4;%if ^^ is yes. if you take out stim shank ==3. if all shanks ==4
downsampleYN=1; %do you want to downsample 1=yes
singleCurrent=6; %plot 6uA results
ploty=1; %Do the plotting
if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end
vec = [100;80;50;30;15;0];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
centroidpos_all=[];
peak_all=[];
samples_rand=[];
singleelectspread=[];
dualelectspread=[];
trials=[1 2 3 4 5];
clear spread_electno spread_electnomaxall sepdistsigdwnsampe sepdistsigdwnsampe_dual
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    dualelectspreadsepdist.(sepcheck)=[];

    hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
    raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
    map = interp1(vec,raw,linspace(100,0,N),'pchip');
    cmap=colormap(map);
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        if length(AMP)>1 && ploty==1
            figure(1+sepdist)
            subplot(2,4,current)
            hold on
        end
        savdwnsamplespkrate=zeros(32,5);
        if ploty==1 && length(AMP)==1
            figure(1+sepdist)
            axes('Position',[0.13         0.112396822842341                     0.775         0.62])
            hold on
            set(gca,'FontSize',12)
        end
        for trial=trials
            check=['T' num2str(trial)];
            if (plotAllShanks==1)
                DataInput2Plot=saveCplit.(sepcheck).(currcheck).(check);% single shank - saveshanksepdist.(sepcheck).(currcheck).(check).D1; all shanks - saveCplit.(sepcheck)
                lengthneeded=size(saveCplit.sep7.C1.T1,2)./numshanksToAvg;% find the smallest number of pairs for any current - we need to down-sample the others to match
                
            else
                DataInput2Plot=saveshanksepdist.(sepcheck).(currcheck).(check).D1;
                lengthneeded=size(saveshanksepdist.(sepcheck).C1.T1.D1,2);% find the smallest number of pairs for any current - we need to down-sample the others to match
            end
            dat=DataInput2Plot;
            dat(isinf(dat))=nan;
            if plotAllShanks==1
                clear pairavgcur erpairavgcur
                for paircount=1:size(dat,2)/numshanksToAvg
                    pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*numshanksToAvg+1:numshanksToAvg* paircount),2);
                end
            else
                pairavgcur=dat;
            end
            
            if downsampleYN==1
                if trial==1
                    size(DataInput2Plot,2)
                    [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
                    samples_rand.(sepcheck)=samples_rand1;
                end
                pairavg_dwnsample=pairavgcur(:, samples_rand.(sepcheck));% downsample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                pairavg_dwnsample=pairavgcur;
            end
            % pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan; %remove rows with less than 10 pairs in the average
            pairavgdwnsample_all_nomean.(sepcheck).(currcheck).(check)=pairavg_dwnsample;
            pairavgdwnsample_all.(sepcheck).(currcheck)(:,trial)=nanmean(pairavg_dwnsample,2);
%             pairavgdwnsample_all.(sepcheck).(currcheck)(:,trial)=nanmean(dat,2);
%             pairavg_dwnsample=dat;

            
            if ploty==1
                stdshade(pairavg_dwnsample',0.2,cmap(trial*floor((length(cmap))/5),:));
            end
            check=['T' num2str(trial)];
         end

        
        if ploty==1
            ylim([0, round(max(pairavgdwnsample_all.(sepcheck).(currcheck),[],'all')/50)*50+50])
            %ylim([0 700])
            if sepdist==5
                 xlim([16-3 16+sepdist+1+3])
            elseif sepdist==7
                 xlim([16-2 16+sepdist+1+2])
            elseif sepdist==9
                 xlim([16-1 16+sepdist+1+1])
            end
            xline(16,'r')
            xline(16+sepdist+1,'r')
            xt = xticks;
            xtl=(xt-16)*50;
            xticklabels(xtl)
            yt = yticks;
            yticklabels(yt(1:end-1));
            ylabel('Firing rate (Sp/s)')
            xlabel('Distance from deepest stim elect (\mum)')
            set(gca,'TickDir','out');
            ax1=gca;
            legend('0:100','25:75','50:50','75:25','100:0')
            title([num2str(AMP(current)) '\muA'])
        end
        
    end
    if ploty==1
        set(gca,'TickDir','out');
        hex = ['#1b0c36';'#532e9e';'#147df5';'#02a612';'#fffb7d';'#ffe26e'];
        raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
        map = interp1(vec,raw,linspace(100,0,N),'pchip');
        cmap=colormap(map);
    end

    
%     if (exist('sepdistsig'))==1 %checks spread and significance of single vs dual
%         dualspread=[];
%         singlespread=[];
%         figure (30)
%         subplot(1,3,round((sepdist-4)/2))
%         hold on
%         for paircount=1:size(sepdistsig.(sepcheck).C6,2)/2
%             sepdistsigdwnsampe.(sepcheck)(:,paircount)=nanmean(nanmean(sepdistsig.(sepcheck).C6(:,( paircount-1)*2+1:2* paircount),2));
%         end
%         for paircount=1:size(sepdistsig_dual.(sepcheck).C6,2)/3
%            sepdistsigdwnsampe_dual.(sepcheck)(:,paircount)=nanmean(nanmean(sepdistsig_dual.(sepcheck).C6(:,( paircount-1)*3+1:3* paircount),2));
%         end
%         sepdistsigdwnsampe_dual.(sepcheck)=[sepdistsigdwnsampe_dual.(sepcheck) nan(1, size(sepdistsigdwnsampe.(sepcheck),2)-size(sepdistsigdwnsampe_dual.(sepcheck),2))];
%         psepsig.(sepcheck)=signrank(sepdistsigdwnsampe_dual.(sepcheck),sepdistsigdwnsampe.(sepcheck),'tail','right');
%         ratiovar=sepdistsigdwnsampe_dual.(sepcheck)./sepdistsigdwnsampe.(sepcheck);
%         ratiovar(isinf(ratiovar))=nan;
%         avgratio.(sepcheck)=nanmean(ratiovar,'all')*100-100;
%     end
%     

    % centroid per pair %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcur=[];
    
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for trial=trials
            check=['T' num2str(trial)];
            dat=pairavgdwnsample_all_nomean.(sepcheck).(currcheck).(check);
            dat(isinf(dat))=nan;
            %dat(sum(~isnan(dat),2)<3,:)=nan;
            pairavgcur=zeros(size(dat,1),size(dat,2));
            clear pairavgcur erpairavgcur
            for paircount=1:size(dat,2)
                pairavgcur(:,paircount)=dat(:,paircount);
                [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur(:,paircount),sepdist);
                centroidpos_all.(sepcheck).(currcheck)(trial,paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%
            end
%             if (plotAllShanks==1)
%                 for paircount=1:size(dat,2)/numshanksToAvg
%                     pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*numshanksToAvg+1:numshanksToAvg* paircount),2);
%                     [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur(:,paircount),sepdist);
%                     centroidpos_all.(sepcheck).(currcheck)(trial,paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%
%                 end
%             else
%                 for paircount=1:size(dat,2)
%                     pairavgcur(:,paircount)=dat(:,paircount);
%                     [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur(:,paircount),sepdist);
%                     centroidpos_all.(sepcheck).(currcheck)(trial,paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%
%                 end
%             end
        end

         %centroidpos_all.(sepcheck).(currcheck)(sum(~isnan(centroidpos_all.(sepcheck).(currcheck)),2)<10,:)=nan;
%         avgpcurr=mean(centroidpos_all.(sepcheck).(currcheck),2);
%         stdavg=std(centroidpos_all.(sepcheck).(currcheck),[],2)./sqrt(size(centroidpos_all.(sepcheck).(currcheck),2));
        %peak_all.(sepcheck).(currcheck)=nan(5,1000);
           %%%%%%%%%%%%%%%%%smoothed peaks
           for trial=trials
               check=['T' num2str(trial)];
               dat=pairavgdwnsample_all_nomean.(sepcheck).(currcheck).(check);
               dat(isinf(dat))=nan;
               %dat(sum(~isnan(dat),2)<3,:)=nan;
               clear pairavgcur erpairavgcur
               for paircount=1:size(dat,2)
                   datnonan=dat(:,paircount);
                   firstnonan=find(~isnan(datnonan),1,'first');
                   lastnonan=find(~isnan(datnonan),1,'last');
                   datnonan(isnan(datnonan))=[];
                   FilterLength=5;
                   b=ones(FilterLength,1)./FilterLength;
                   smoothedenvelope=filtfilt(b,1,datnonan);
                   %smoothedenvelope([1 end],1)=nan;
                   %sm1 = conv(datnonan,b,'same');
                   pairavgcur(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];
                   %                    peak=find(pairavgcur(:,paircount)==max(pairavgcur(:,paircount),[],'omitnan'));
                   %                    %[electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur(:,paircount),sepdist);
                   %                    peak_all.(sepcheck).(currcheck)(trial,paircount)=(peak(1)-16).*50;
                   %plot(pairavgcur(:,paircount))
               end
               %figure(trial+sepdist*100); hold on; 
%                if (plotAllShanks==1)
%                    for paircount=1:size(dat,2)/numshanksToAvg
%                        datnonan=nanmean(dat(:,( paircount-1)*numshanksToAvg+1:numshanksToAvg* paircount),2);
%                        firstnonan=find(~isnan(datnonan),1,'first');
%                        lastnonan=find(~isnan(datnonan),1,'last');
%                        datnonan(isnan(datnonan))=[];
%                        FilterLength=5;
%                        b=ones(FilterLength,1)./FilterLength;
%                        smoothedenvelope=filtfilt(b,1,datnonan);
%                        pairavgcur(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];
%                    end
%                else
%                    for paircount=1:size(dat,2)
%                        datnonan=dat(:,paircount);
%                        firstnonan=find(~isnan(datnonan),1,'first');
%                        lastnonan=find(~isnan(datnonan),1,'last');
%                        datnonan(isnan(datnonan))=[];
%                        FilterLength=5;
%                        b=ones(FilterLength,1)./FilterLength;
%                        smoothedenvelope=filtfilt(b,1,datnonan);
%                        %smoothedenvelope([1 end],1)=nan;
%                        %sm1 = conv(datnonan,b,'same');
%                        pairavgcur(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];
%                        %                    peak=find(pairavgcur(:,paircount)==max(pairavgcur(:,paircount),[],'omitnan'));
%                        %                    %[electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur(:,paircount),sepdist);
%                        %                    peak_all.(sepcheck).(currcheck)(trial,paircount)=(peak(1)-16).*50;
%                        %plot(pairavgcur(:,paircount))
%                    end
%                end
               %pairavgcur(sum(~isnan(pairavgcur),2)<25,:)=nan;
               pairavgcur(:,sum(pairavgcur>0,1)==0)=nan;
               [peak,c]=find(pairavgcur==max(pairavgcur));
               %peak=find(pairavgcur(:,paircount)==max(pairavgcur(:,paircount),[],'omitnan'));
               peak_all.(sepcheck).(currcheck)(trial,c)=(peak'-16).*50;
%                lengthneeded=size(saveshanksepdist.(sepcheck).C1.T1.D1,2);% find the smallest number of pairs for any current - we need to down-sample the others to match
%                if trial==1
%                    [samples_rand1]=DownSample(dat,lengthneeded,s,seedpoint);
%                    samples_rand.(sepcheck)=samples_rand1;
%                end
%                figure(sepdist+50)
%                hold on
%                stdshade(pairavgcur',0.2,cmap(trial*floor((length(cmap))/5),:));
           end
           if ploty==1
               if sepdist==5
                   xlim([16-3 16+sepdist+1+3])
               elseif sepdist==7
                   xlim([16-2 16+sepdist+1+2])
               elseif sepdist==9
                   xlim([16-1 16+sepdist+1+1])
               end
               xline(16,'r')
               xline(16+sepdist+1,'r')
               xt = xticks;
               xtl=(xt-16)*50;
               xticklabels(xtl)
               ylabel('Firing rate (Sp/s)')
               xlabel('Distance from deepest stim elect (\mum)')
               set(gca,'TickDir','out');
               legend('0:100','25:75','50:50','75:25','100:0')
           end
%            peak_all.(sepcheck).(currcheck)=peak_all.(sepcheck).(currcheck)(:, samples_rand.(sepcheck));
%            peak_all.(sepcheck).(currcheck)(sum(~isnan(peak_all.(sepcheck).(currcheck)),2)<10,:)=nan;
%     
           %%%%%%%%%%%%%%%%%%%%plotting
           avgpcurrc=mean(centroidpos_all.(sepcheck).(currcheck),2);
        stdavgc=std(centroidpos_all.(sepcheck).(currcheck),[],2)./sqrt(size(centroidpos_all.(sepcheck).(currcheck),2));
            avgpcurr=nanmean(peak_all.(sepcheck).(currcheck),2);
            %avgpcurr=peak_all.(sepcheck).(currcheck)(:,col2plot);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            stdavg=nanstd(peak_all.(sepcheck).(currcheck),[],2)./sqrt(sum(~isnan(peak_all.(sepcheck).(currcheck)),2));
            %stdavg=0;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if ploty==1
            avgpcurr=avgpcurr(trials);
            if length(AMP)==1
            figure(1+sepdist)
            axes('Position',[ax1.Position(1) .73 ax1.Position(3) .2])
            ax2=gca;
            hold on
            set(gca,'FontSize',12)
            for colorbar=1:length(avgpcurr)
                color_current=cmap(trials(colorbar)*floor((length(cmap))/5),:);
                er = errorbar(avgpcurr(colorbar),colorbar,stdavg(trials(colorbar)),stdavg(trials(colorbar)),'horizontal');er.Color = [0, 0, 0]; er.LineStyle = 'none';
                scatter(avgpcurr(colorbar), colorbar, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)
                %er = errorbar(avgpcurrc(colorbar),colorbar+5,stdavgc(trials(colorbar)),stdavgc(trials(colorbar)),'horizontal');er.Color = [0, 0, 0]; er.LineStyle = 'none';
                %scatter(avgpcurrc(colorbar), colorbar+5, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)
            end
            yticks([0 1 2 3 4 5 6 7 8])
            labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
            labels=labels(trials);
            yticklabels([{''} labels {''}])
            ylim([0.5 length(trials)+0.5])
            xlim((ax1.XLim-16).*50)
            set(gca,'TickDir','out');
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            xline(0,'r')
            xline((sepdist+1)*50,'r')
            end
            
            
         

            
        end
    end
end



%% Model
%C:\Users\smei0006\Documents\Experimental_Design\predictDualEnvelope.m
%Data input can be found E:\DATA\savecsplit.mat 
%Alternatively just run it all again up above but will take a long time

%5050%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Additive not scaled - make sure you have half amplitude data!
%saveCplit=saveCplit_all4shanks_half;
%Hypoth 1: multiply1=1; 

%Weighted not scaled - make sure you have full amplitude data!
%saveCplit=saveCplit_all4shanks;
%Hypoth 1: multiply1=0.5; 

%Additive scaled - make sure you have half amplitude data!
%saveCplit=saveCplit_all4shanks_half;
%Hypoth 4/5

%Weighted scaled - make sure you have full amplitude data!
%saveCplit=saveCplit_all4shanks;
%Hypoth 4/5



%7525%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%got to 7525 section line 640

%Additive not scaled - make sure you have quarter amplitude data!
%saveCplit=saveCplit_all4shanks_quarter;
%Hypoth 1: which75_25condition=2; 

%Weighted not scaled - make sure you have full amplitude data!
%saveCplit=saveCplit_all4shanks_quarter;
%Hypoth 2: which75_25condition=2; for plotting and 1 will need to be run for centroid error 

%Additive scaled - make sure you have quarter amplitude data!
%saveCplit=saveCplit_all4shanks_quarter;
%Hypoth 6/7/10: which75_25condition=2;

%Weighted scaled - make sure you have full amplitude data!
%saveCplit=saveCplit_all4shanks;
%Hypoth 6/7/10: which75_25condition=2; for plotting and 1 will need to be run for centroid error 


%all combined to get two sclaing factors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%saveCplit=saveCplit_all4shanks;
%altogether model #2: mdlcoef are the coefficients

%% model shift as a percentage of sep dist
mean1=0;
mean2=0;
colorchoice={'r' 'b' 'k'};
for sepdist=5:2:9
    checksep=['sep' num2str(sepdist)];
b1centroid.(checksep)=[b1centroid_2575.(checksep); b1centroid_5050.(checksep); b1centroid_7525.(checksep);];
b2centroid.(checksep)=[b2centroid_2575.(checksep); b2centroid_5050.(checksep); b2centroid_7525.(checksep);];

figure (1)
hold on
diffcentroid=abs(flipud(centroidpos_all.(checksep).C6(2:4,:))-b1centroid.(checksep));
amean=nanmean(diffcentroid,2);
astd=nanstd(diffcentroid,[],2)./sqrt(sum(~isnan(diffcentroid),2)); 
er = errorbar(2:4,amean,astd);er.Color = colorchoice{sepdist/2-1.5}; er.LineStyle = 'none';
scatter(2:4,amean,[], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', colorchoice{sepdist/2-1.5})
mean1=mean1+mean(diffcentroid,'all');

figure (2)
hold on
diffcentroid=abs(flipud(centroidpos_all.(checksep).C6(2:4,:))-b2centroid.(checksep));
amean2=nanmean(diffcentroid,2);
astd2=nanstd(diffcentroid,[],2)./sqrt(sum(~isnan(diffcentroid),2)); 
er = errorbar(2:4,amean,astd);er.Color = colorchoice{sepdist/2-1.5}; er.LineStyle = 'none';
scatter(2:4,amean,[], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', colorchoice{sepdist/2-1.5})
mean2=mean2+mean(diffcentroid,'all');
end
mean1/3
mean2/3
figure(1)
xticks(1:5)
labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
labels=labels(trials);
xticklabels(labels)
ylabel('Peak prediction model error (\mum)')
xlabel('Percentage current delivered')
ylim([0 75])
xlim([1 5])
set(gca,'TickDir','out');
legend('300\mum', '400\mum', '500\mum')
title('B1')

figure (2)
xticks(1:5)
labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
labels=labels(trials);
xticklabels(labels)
ylabel('Peak prediction model error (\mum)')
xlabel('Percentage current delivered')
ylim([0 75])
xlim([1 5])
set(gca,'TickDir','out');
legend('300\mum', '400\mum', '500\mum')
title('B2')



%%  additional significance tests

for current=2:length(AMP)
    %is the spike rate significantly different for single vs dual
    currcheck=['C' num2str(AMP(current))];
    sep5spkdual=(max(pairavgdwnsample_all_nomean.sep5.(currcheck).T2)+max(pairavgdwnsample_all_nomean.sep5.(currcheck).T3)+max(pairavgdwnsample_all_nomean.sep5.(currcheck).T4))./3;
    sep9spkdual=(max(pairavgdwnsample_all_nomean.sep9.(currcheck).T2)+max(pairavgdwnsample_all_nomean.sep9.(currcheck).T3)+max(pairavgdwnsample_all_nomean.sep9.(currcheck).T4))./3;
    sep7spkdual=(max(pairavgdwnsample_all_nomean.sep7.(currcheck).T2)+max(pairavgdwnsample_all_nomean.sep7.(currcheck).T3)+max(pairavgdwnsample_all_nomean.sep7.(currcheck).T4))./3;
    sep5spksingle=(max(pairavgdwnsample_all_nomean.sep5.(currcheck).T1)+max(pairavgdwnsample_all_nomean.sep5.(currcheck).T5))./2;
    sep9spksingle=(max(pairavgdwnsample_all_nomean.sep9.(currcheck).T1)+max(pairavgdwnsample_all_nomean.sep9.(currcheck).T5))./2;
    sep7spksingle=(max(pairavgdwnsample_all_nomean.sep7.(currcheck).T1)+max(pairavgdwnsample_all_nomean.sep7.(currcheck).T5))./2;
    
    [p55_spk,h,stats] = signrank(sep5spkdual(:), sep5spksingle(:),'tail','right');
    [p77_spk,h,stats] = signrank(sep7spkdual(:), sep7spksingle(:),'tail','right');
    [p99_spk,h,stats] = signrank(sep9spkdual(:), sep9spksingle(:),'tail','right');
    p55_spk_all(current-1)=p55_spk;
    p77_spk_all(current-1)=p77_spk;
    p99_spk_all(current-1)=p99_spk;
end

% is the centroid shift significant 
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        [p,tbl,stats]=anova2(centroidpos_all.(sepcheck).(currcheck)',1,'off'); % for parametric
        pall.(sepcheck) (current)=p(1);
        [p,tbl,stats]=friedman(centroidpos_all.(sepcheck).(currcheck)',1,'off'); %for non-para
        pall2.(sepcheck) (current)=p;
    end
end

%%
% is the peak shift significant 
for sepdist=9%5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        [p,tbl,stats]=anova2(peak_all.(sepcheck).(currcheck)',1,'on'); % for parametric
        pall.(sepcheck) (current)=p(1);
        [p,tbl,stats]=friedman(peak_all.(sepcheck).(currcheck)',1,'off'); %for non-para
        [c,~,~,gnames] = multcompare(stats,'Alpha',0.01);
        c(:,7)=zeros(size(c,1),1);
        for i=1:size(c,1)
            c(i,7)=ranksum(peak_all.(sepcheck).(currcheck)(c(i,1),:)', peak_all.(sepcheck).(currcheck)(c(i,2),:)');
        end
        pall2.(sepcheck) (current)=p;
    end
end