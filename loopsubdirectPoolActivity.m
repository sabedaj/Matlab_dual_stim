marmoflag=0;

startpointseconds=2;
secondstoanalyse=8;
sepdist=-6;
for sepdist=-6:-2:-10
%%% Laminar
overallpoolact=[];
meananimal=[];
overallcentroidshift=zeros(5,2000);
numrejectedchannels=0;
numacceptchan=0;
AMP=[0 1 2 3 4 6 8 10];

for trial=1:5
    check1=['T' num2str(trial)];
    overallpoolact.(check1)=[];
    overallpool_nostim.(check1)=[];
    for current = 1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        Csplit.(currcheck).(check1)=[];
        spread.(currcheck).(check1)=[];
        spkrateatcentroid.(currcheck).(check1)=[];
        CsplitCentroid.(currcheck)=[];
         for shanksep=1:3
             bincheck=['D' num2str(shanksep)];
            bin.(currcheck).(check1).(bincheck)=[];
         end
    end
end
counternumberofpairs=0;
countershank=0;
stimshank=[];
AR=[1 2 3 4];%array shanks
binchecksave='blank';

for ratN=14:20
if ratN<10 && marmoflag==0
    Ratnum=['Rat_00' num2str(ratN)];
elseif ratN>10 && marmoflag==0
    Ratnum=['Rat_0' num2str(ratN)];
else%%marmo
    Ratnum=['Marmo_00' num2str(ratN)];
end
%%%% Asplit for animal averaging
AMP=[0 1 2 3 4 6 8 10 15];
for trial=1:5
    check1=['T' num2str(trial)];
    for current = 1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        Asplit.(Ratnum).(currcheck).(check1)=[];
         AsplitCentroid.(Ratnum).(currcheck)=[];
    end
end



%%%%
cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
D_data=dir;

savestimshanks=zeros(length(D_data),1);
for k = 3:length(D_data) % avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
%         if strcmp(D_data(k).name,('S4E4_9elect_001_210325_134828'))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             continue
%         end
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


%         %% timing testing
%         trig = loadTrig(0);
%         Overall_time_to_analyse=0;%time from beginning of Startpoint_analyse (remembering there is 20s of no stim at beginning) %set to zero fo no input
%         startpointseconds=8; %How long after the trigger do you want skip spike analysis(ms)?
%         secondstoanalyse=10; %How long after the trigger do you want to analyse spikes for(ms)?
%         starttrial=1;
%         trialjump=1;
%         TrialParams=loadTrialParams;
%         endtrial=max(cell2mat(TrialParams(:,2)));
%         printspiking=0;
%         [IDstruct, baslinespikestruct]=sortTrials_SM(startpointseconds,secondstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,Overall_time_to_analyse);
%         %             load('IDstruct.mat','IDstruct')%%%%%%% testing no basline removal
%         [avgnospT,~,~] = AverageTrialResponse_SM(IDstruct,baslinespikestruct);%%%%% testing no basline removal or timing
        
        [normspk_all,centroidpershank]=PoolNormalisedActivity(avgnospT,startpointseconds, secondstoanalyse);
        if sum(normspk_all.S1_T1,'all')==0 && sum(centroidpershank,'all')==140
            continue
        end
        if stimChn(1)<17
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
        AMP=loadAMP;
        AMP(AMP==-1)=0;
        for shanknum=1:4
            if shanknum==shank
                savestimshanks(k)=shanknum;
                continue
            end
            
%             if shanknum~=shank+1 || shanknum~=shank-1
% 
%                 continue
%             end
            countershank=countershank+1;
            trialcheck=0;
            overallcentroidshift(1:5,(counternumberofpairs)*3*length(AMP)+(countershank-1)*(length(AMP)-1)+countershank:(counternumberofpairs)*3*length(AMP)+(countershank*length(AMP)))=centroidpershank((shanknum-1)*5+1:(shanknum*5),:);
            for trial=1:5
                check=['S' num2str(shanknum) '_T' num2str(trial)];
                rejectchan=find(normspk_all.(check)==0);
                numrejectedchannels=length(rejectchan)+numrejectedchannels;
                acceptchan=find(normspk_all.(check)~=0);
                numacceptchan=length(acceptchan)+numacceptchan;
                check1=['T' num2str(trial)];
                overallpoolact.(check1)=[overallpoolact.(check1),[zeros(16-stimChn(1),size(normspk_all.(check),2)-1); normspk_all.(check)(:,2:end); zeros(32-(16-stimChn(1))-16,size(normspk_all.(check),2)-1)]];
                for current = 1:length(AMP)
                    currcheck=['C' num2str(AMP(current))];
                    Csplit.(currcheck).(check1)=[Csplit.(currcheck).(check1) [zeros(16-stimChn(1),1); normspk_all.(check)(:,current); zeros(32-(16-stimChn(1))-16,1)]];
                    spread.(currcheck).(check1)=[spread.(currcheck).(check1); numel(find(normspk_all.(check)(:,current)>(sqrt(var(normspk_all.(check)(:,current)))+mean(normspk_all.(check)(:,current)))))];
                    if centroidpershank((shanknum-1)*5+trial,current)-1<=0
                        spkrateatcentroid.(currcheck).(check1)=[spkrateatcentroid.(currcheck).(check1); mean(normspk_all.(check)(centroidpershank((shanknum-1)*5+trial,current):centroidpershank((shanknum-1)*5+trial,current)+1,current))];
                    else
                        spkrateatcentroid.(currcheck).(check1)=[spkrateatcentroid.(currcheck).(check1); mean(normspk_all.(check)(centroidpershank((shanknum-1)*5+trial,current)-1:centroidpershank((shanknum-1)*5+trial,current)+1,current))];
                    end
                    Asplit.(Ratnum).(currcheck).(check1)=[Asplit.(Ratnum).(currcheck).(check1) [zeros(16-stimChn(1),1); normspk_all.(check)(:,current); zeros(32-(16-stimChn(1))-16,1)]];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if trialcheck~=1
                        CsplitCentroid.(currcheck)=[CsplitCentroid.(currcheck) centroidpershank((shanknum-1)*5+1:(shanknum*5),current)-stimChn(1)];
                        AsplitCentroid.(Ratnum).(currcheck)=[AsplitCentroid.(Ratnum).(currcheck) centroidpershank((shanknum-1)*5+1:(shanknum*5),current)];%%%%%%%%%%%%%%%%%%%%
                        if isfield(meananimal, Ratnum) && isfield(meananimal.(Ratnum), currcheck) && countershank==1
                            meananimal.(Ratnum).(currcheck)=meananimal.(Ratnum).(currcheck)+1;
                        elseif countershank==1
                            meananimal.(Ratnum).(currcheck)=1;
                        end
                        
                    end

                   % plot resutls of 200um away, 400um away and 600um away
                        Shanksep=shanknum-shank;
                        bincheck=['D' num2str(abs(Shanksep))];
                        if strcmp(bincheck,binchecksave)
                            bin.(currcheck).(check1).(bincheck)(:,end)=mean([bin.(currcheck).(check1).(bincheck)(:,end) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]],2);
                        else
                            bin.(currcheck).(check1).(bincheck)=[bin.(currcheck).(check1).(bincheck) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                        end

                end
                overallpool_nostim.(check1)=[overallpool_nostim.(check1),[zeros(16-stimChn(1),1); normspk_all.(check)(:,1); zeros(32-(16-stimChn(1))-16,1)]];
                trialcheck=1;

            end
            binchecksave=bincheck;
        end
        counternumberofpairs=counternumberofpairs+1;
        countershank=0;
        binchecksave='blank';
%         for current=1:length(AMP)
%             currcheck=['C' num2str(AMP(current))];
%             if k~=3
%                 currcheck_prev=['C' num2str(AMP(current-1))];
%                 meananimal.(Ratnum).(currcheck)=size(CsplitCentroid.(currcheck),2)-meananimal.(Ratnum).(currcheck_prev);
%             else
%                 meananimal.(Ratnum).(currcheck)=size(CsplitCentroid.(currcheck),2);
%             end
%         end
end
savestimshanks(savestimshanks==0)=[];
stimshank=[stimshank; savestimshanks];
end
overallcentroidshift_nostim=overallcentroidshift(:,1:8:counternumberofpairs*3*8);
overallcentroidshift_stimonly=overallcentroidshift;
overallcentroidshift_stimonly(:,1:8:counternumberofpairs*3*8)=[];
data=diff(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7));
check=['sep' num2str(-1*sepdist-1)];
saveCplit.(check)=Csplit;
saveCsplitcentroid.(check)=CsplitCentroid;
savespkrateatcentroid.(check)=spkrateatcentroid;
savespread.(check)=spread;
end
%%
% AR=[1 2 3 4];%array shanks
% AMP=[0 1 2 3 4 6 8 10];
% clear pairavg stdpairavg cavgspk
% for loopamp=1:length(AMP)
%     currcheck=['C' num2str(AMP(loopamp))];
%     for looppair=1:length(stimshank)
%     Shanksep=AR-stimshank(looppair);
%     Shanksep(Shanksep==0)=[];
%     bincheck=['D' num2str(abs(Shanksep(1)))];
%     bin.(currcheck).(bincheck)=[bin.(currcheck).(bincheck) CsplitCentroid.(currcheck)(:,(looppair-1)*3+1)];
%         bincheck=['D' num2str(abs(Shanksep(2)))];
%     bin.(currcheck).(bincheck)=[bin.(currcheck).(bincheck) CsplitCentroid.(currcheck)(:,(looppair-1)*3+2)];
%         bincheck=['D' num2str(abs(Shanksep(3)))];
%     bin.(currcheck).(bincheck)=[bin.(currcheck).(bincheck) CsplitCentroid.(currcheck)(:,(looppair-1)*3+3)];
%     
%     end
%     
% end
%%

AMP=[0 1 2 3 4 6 8 10];
clear pairavg stdpairavg cavgspk 
for looppair=1:length(AMP)
    currcheck=['C' num2str(AMP(looppair))];
    
% Define number of columns to average
AVG_COLS = 3;
% Dimension over which to average
DIM = 2; % Columns
% Use filter to calculate the moving average across EVERY combination of columns
pairavg.(currcheck) = filter(ones(1,AVG_COLS)/AVG_COLS,1,CsplitCentroid.(currcheck).*50,[],DIM);
% Grab only the column averages that were actually wanted
pairavg.(currcheck) = pairavg.(currcheck)(:,AVG_COLS:AVG_COLS:end);
% 
% for tsplit=1:5
%     T=['T' num2str(tsplit)];
%     cplitavg.(currcheck).(T) = filter(ones(1,AVG_COLS)/AVG_COLS,1,Csplit.(currcheck).(T),[],DIM);
%     cplitavg.(currcheck).(T) =  cplitavg.(currcheck).(T)(:,AVG_COLS:AVG_COLS:end);
%     cplitavg.(currcheck).(T)(cplitavg.(currcheck).(T)==0)=nan;
%     if tsplit==1
%         cavgspk.(currcheck)=nanmean(cplitavg.(currcheck).(T),2);
%     else
%         cavgspk.(currcheck)=[cavgspk.(currcheck) nanmean(cplitavg.(currcheck).(T),2)];
%     end
% end


%stdpairavg(:,looppair)=std(CsplitCentroid.(currcheck).*50,[],2)./sqrt(3);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%should I be calculating this here?
%counternumberofpairs*3*length(AMP)).*50,2);
%stpairacg(looppair,:)=std(overallcentroidshift_stimonly(:,(looppair-1)*3+1:3:counternumberofpairs*3*length(AMP)).*50,2)./sqrt(
end
%%
for asplit=1:length(AMP)
    csp=['C' num2str(AMP(asplit))];
currentsplit=pairavg.(csp)';
[p,tbl,stats] = kruskalwallis(currentsplit);
cp(asplit,1)=p;
end
%%
figure
for current=1:length(AMP)
subplot(2,4,current)
    hold on
    currcheck=['C' num2str(AMP(current))];
    avgpcurr=mean(pairavg.(currcheck),2);%mean(pairavg(:,current:length(AMP):counternumberofpairs*length(AMP)),2);
    stdavg=std(pairavg.(currcheck),[],2)./sqrt(size(pairavg.(currcheck),2));

    bar(avgpcurr)
er = errorbar(1:5,avgpcurr,stdavg,stdavg);er.Color = [0 0 0]; er.LineStyle = 'none';
%er2(current,1:5)=std(overallcentroidshift_stimonly(:,current:length(AMP):counternumberofpairs*3*length(AMP)).*50,[],2)./sqrt(counternumberofpairs*3);
ylabel('Distance from deepest stim electrode (\mum)')
xlabel('Percentage current delivered')
xticks([0 1 2 3 4 5 6 7 8])
xticklabels(fliplr({'' '0:100','25:75','50:50','75:25','100:0' ''}))
title(['Centroid location: ' num2str(AMP(current)) '\muA total'])
xlim([0.5 5.5])
ylim([-50 400])
end

%% figure for split shanks
AMP=[0 1 2 3 4 6 8 10 30];
cmap=colormap('gray');
linS = {'.','o','+','x','*'};

for current=1:length(AMP)
    figure
    currchec=['C' num2str(AMP(current))];
    
    
    for shankit=1:3
        dcheck=['D' num2str(shankit)];
        subplot(1,3,shankit)
        for trial=1:5
            hold on
            check=['T' num2str(trial)];
            
            dat=bin.(currchec).(check).(dcheck);
            %dat(dat==0)=nan;
            dat(sum(~isnan(dat),2)<1,:)=nan;
            pairavgcur=nanmean(dat,2);
            erpairavgcur=nanstd(dat,[],2);
            
            %find(~isnan(nanmean(pairavgcur,2)),'first');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plot(nanmean(pairavgcur,2),'Color',cmap(trial*floor((length(cmap)-50)/5),:))
            errorbar(nanmean(pairavgcur,2),nanmean(erpairavgcur,2),'Color',cmap(trial*floor((length(cmap)-50)/5),:),'marker', linS{trial})
            % dat(sum(~isnan(dat),2)<6,:)=nan;
            stder=nanstd(dat,[],2)./sqrt(sum(~isnan(dat),2));
            s(trial,current)=nanmean(stder);
            %errorbar(nanmean(dat,2),stder,'Color',cmap(trial*floor((length(cmap)-50)/5),:),'marker', linS{trial})
        end
            title([num2str(AMP(current)) '\muA: Shank distance ' num2str(shankit*200) '\mum'])
    ylim([0, 600])
    xline(16,'r')
    xline(16-sepdist,'r')
    xt = xticks;
    xtl=(xt-16+(sepdist/2))*50;
    xticklabels(xtl)
    ylabel('Sp/s')%'Normalised spike rate'
    xlabel('Distance from center of electrode pair (\mum)')
    xlim([15 27])
    end
   
end

%%
%Csplit=saveCplit.sep5;
spreadv2=zeros(length(AMP),5);
AMP=[0 1 2 3 4 6 8 10];
cmap=colormap('gray');
linS = {'.','o','+','x','*'};
figure
for current=1:length(AMP)
    currcheck=['C' num2str(AMP(current))];
    subplot(2,5,current)
    for trial=1:5
        hold on
        check=['T' num2str(trial)];
        dat=Csplit.(currcheck).(check);
        dat(dat==0)=nan;
        dat(sum(~isnan(dat),2)<5,:)=nan;
        clear pairavgcur erpairavgcur
        for paircount=1:size(dat,2)/3
            pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*3+1:3* paircount),2);
           % erpairavgcur(:,paircount)=nanstd(dat(:,( paircount-1)*3+1:3* paircount),[],2)./sqrt(sum(~isnan(dat(:,( paircount-1)*3+1:3* paircount)), 2));
        end
        %find(~isnan(nanmean(pairavgcur,2)),'first');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plot(nanmean(pairavgcur,2),'Color',cmap(trial*floor((length(cmap)-50)/5),:))
        
        avgpair=nanmean(pairavgcur,2);
       spreadv2(current,trial)=numel(find(avgpair>max(avgpair)*0.5));
        %errorbar(nanmean(pairavgcur,2),nanmean(erpairavgcur,2),'Color',cmap(trial*floor((length(cmap)-50)/5),:),'marker', linS{trial})
       % dat(sum(~isnan(dat),2)<6,:)=nan;
        stder=nanstd(dat,[],2)./sqrt(sum(~isnan(dat),2));
        s(trial,current)=nanmean(stder);
        %errorbar(nanmean(dat,2),stder,'Color',cmap(trial*floor((length(cmap)-50)/5),:),'marker', linS{trial})
        
    end
    title([num2str(AMP(current)) '\muA'])
    ylim([0, 200])
    xline(16,'r')
    xline(16-sepdist,'r')
    xt = xticks;
    xtl=(xt-16+(sepdist/2))*50;
    xticklabels(xtl)
    ylabel('Sp/s')
    xlabel('Relative electrode number')
end

%%
cmap=colormap('gray');
figure
hold on
for current=2:length(AMP)
plot(spreadv2(current,1:2:end),'Color',cmap(current*floor((length(cmap)-50)/length(AMP)),:))
end
%%
%shuffle
shufperm=reshape(randperm(5*counternumberofpairs*3*7), 5, counternumberofpairs*3*7);
suffle_permutation=overallcentroidshift_stimonly(shufperm);
figure 
hold on
bar(mean(suffle_permutation(:,1:counternumberofpairs*3*7).*50,2))
ylabel('Distance from deepest electrode (\mum)')
xlabel('Percentage current delivered')
xticklabels(fliplr({'','0:100','25:75','50:50','75:25','100:0',''}))
title('Centroid location: shuffled data')
[p,tbl,stats]=anova1((suffle_permutation(:,1:counternumberofpairs*3*7))'.*50);
ylabel('Distance from deepest electrode (\mum)')
xlabel('Percentage current delivered')
xticklabels(fliplr({'0:100','25:75','50:50','75:25','100:0'}))
title('Centroid location: shuffled data')

figure 
hold on
bar(mean(diff(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7)).*-50,2))
er = errorbar(1:4,mean(diff(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7)).*-50,2),std(diff(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7)).*-50,[],2)./sqrt(counternumberofpairs*3*7),...
    std(diff(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7)).*-50,[],2)./sqrt(counternumberofpairs*3*7));er.Color = [0 0 0]; er.LineStyle = 'none';
ylabel('Average shift per 25% change in distribution of current (\mum)')
xlabel('Percentage current delivered')
xticklabels(fliplr({'','0:100','25:75','50:50','75:25','100:0',''}))
title('centroid location')
avgnostim=mean(diff(overallcentroidshift_nostim(:,1:counternumberofpairs*3)),'all')*-1*50;
stder=std(diff(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7)),0,'all')*50/sqrt(length(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7)));
sound(sin(1:3000));
%%
prevfirstelem=1;
prevlastelem=1;
ratePool=[];
rpoolnom=[];
stderall=[];
SMOOTHING=1;
window = [1/4 1/4 1/4 1/4];%normpdf(-2*SMOOTHING:2*SMOOTHING,0,SMOOTHING);
for trial=1:5
    check1=['T' num2str(trial)];
    overallpoolact.(check1)(overallpoolact.(check1)==0)=NaN;
    trialdata=overallpoolact.(check1);

    nonan_len=zeros(size(trialdata,1),1);
    for i=1:size(trialdata,1)
        testlength=trialdata(i,:);
        nonan_len(i) = length(testlength(~isnan(testlength)));
        leng=65; %set to number of trials to avoid edge effects
        nonan_len(nonan_len<leng)=0;
         if nonan_len(i)==0
             trialdata(i,:)=NaN;
            continue
         end
        
%        leng=1; %set to number of trials to avoid edge effects
%         nonan_len(nonan_len<leng)=0;
%          if nonan_len(i)==0
%              trialdata(i,:)=NaN;
%             continue
%          end
%          input=find(~isnan(testlength));
%          s = RandStream('mlfg6331_64','seed','shuffle');
%         [y1,idx] = datasample(s,input,(nonan_len(i)-leng),'Replace',false);
%         trialdata(i,y1)=NaN;
    end
    lengthall = size(overallpoolact.(check1),2);
    if any(diff(find(nonan_len>0))>1)
        testlength=find(nonan_len>0);
         trialdata(testlength((diff(testlength)>1)),:)=NaN;
    end


    trPool=nanmean(trialdata,2);
    trPool(isnan(trPool))=[];
    stder=nanstd(trialdata,0,2)./sqrt(sum((~isnan(trialdata)),2));
    stder(isnan(stder))=[];
        %rate=trPool;
    
    checkfirstelem=find(~isnan(mean(trialdata','omitnan')),1);
    checklastelem=33-find(~isnan(fliplr(mean(trialdata','omitnan'))),1);
    if trial ~=1
        if checkfirstelem>prevfirstelem
            removerow=checkfirstelem-prevfirstelem;
            rpoolnom(1:removerow,:)=[];
            stderall(1:removerow,:)=[];
            prevfirstelem=checkfirstelem;
        elseif checkfirstelem<prevfirstelem
            removerow=prevfirstelem-checkfirstelem;
            trPool(1:removerow,:)=[];
            stder(1:removerow,:)=[];
        else
            prevfirstelem=checkfirstelem;
        end
        if checklastelem>prevlastelem
            removerow=checklastelem-prevlastelem;
            trPool((end-removerow+1):end,:)=[];
            stder((end-removerow+1):end,:)=[];
        elseif checklastelem<prevlastelem
            removerow=prevlastelem-checklastelem;
            rpoolnom((end-removerow+1):end,:)=[];
            stderall((end-removerow+1):end,:)=[];
            prevlastelem=checklastelem;
        else
            prevlastelem=checklastelem;
        end
    else 
        prevlastelem=checklastelem;
        prevfirstelem=checkfirstelem;
    end
    rate=trPool;
    %rate = conv(trPool,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
    %rate =rate(4:end-3);% rate(3*SMOOTHING+1:end-3*SMOOTHING);
    rpoolnom=[rpoolnom,rate];
    %atePool=[ratePool,rate];
    stderall=[stderall,stder];
end


% fig = gcf;
% axObjs = fig.Children;
% dataObjs = axObjs.Children;
% L1=dataObjs(1).YData;
% L2=dataObjs(2).YData;
% L3=dataObjs(3).YData;
% L4=dataObjs(4).YData;
% L5=dataObjs(5).YData;
% 
% E1=dataObjs(3).YNegativeDelta;
% E2=dataObjs(4).YPositiveDelta;
% E3=dataObjs(5).YPositiveDelta;
% E4=dataObjs(6).YPositiveDelta;
% E5=dataObjs(7).YPositiveDelta;
% eall=[E1;E2;E3;E4;E5];
% eall=flipud(eall);
% lall=[L1;L2;L3;L4;L5];
% lall=flipud(lall);
% 
% llinall=[L1;L2;L3;L4;L5];
% llinall=flipud(llinall);

figure
hold on
cmap=colormap('gray');
linS = {'.','o','+','x','*'};
% for i=1:5
% % Fit line to data using polyfit
% c = polyfit(1:size(rpoolnom,1),rpoolnom(:,i),1);
% % Display evaluated equation y = m*x + b
% disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
% % Evaluate fit equation using polyval
% y_est = polyval(c,1:size(rpoolnom,1));
% plot(y_est,'Color',cmap(i*floor((length(cmap)-50)/5),:))
% end
title('Shifting centroid laminar')
xlabel('Relative electrode number')
ylabel('Linear fit response')
xline(16-prevfirstelem+1,':k')
xline(16-prevfirstelem-sepdist+1,':k')

%  figure
% hold on
% windowWidth = 2; 
% window = ones(windowWidth,1) / windowWidth;
% 
% 
% for i=1:5
% rpoolnom_smooth=conv(rpoolnom(:,i),window,'same');
% %rpoolnom_smooth=rpoolnom_smooth(length(window):length(rpoolnom_smooth)-length(window)+1);
% plot(rpoolnom_smooth);
% hold on
% end
% plot(rpoolnom(:,1))
for i=1:5
errorbar(rpoolnom(:,i)',stderall(:,i)','Color',cmap(i*floor((length(cmap)-50)/5),:),'marker', linS{i})%, 'LineStyle','none')
end
% for i=1:5
% plot(1:length(llinall(i,:)),llinall(i,:),'Color',cmap(i*floor((length(cmap)-50)/5),:))
% end
% for i=1:5
%     hold on
% errorbar(lall(i,:),eall(i,:),'Color',cmap(i*floor((length(cmap)-50)/5),:),'marker', linS{i},'LineStyle','none')
% end
title('Shifting centroid laminar')
xlabel('Relative electrode number')
ylabel('Normalised response')
xline(16-prevfirstelem+1,':k')
xline(16-prevfirstelem-sepdist+1,':k')

%%%need to normalise, colour change, show plot with same axes and all
%%%current levels. standard deviation. Only 18 samples since 3 shanks x 6
%%%pairs = 18. very limited maybe reduce to 10uA only - strongest response
%%%approximately at saturation.


%%
figure
hold on
for i=1:5
% Fit line to data using polyfit
c = polyfit(1:size(rpoolnom,1),rpoolnom(:,i),1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
% Evaluate fit equation using polyval
y_est = polyval(c,1:size(rpoolnom,1));
plot(y_est,'Color',cmap(i*floor((length(cmap)-50)/5),:))
end
title('Shifting centroid laminar')
xlabel('Relative electrode number')
ylabel('Linear fit response')
xline(16-prevfirstelem+1,':k')
xline(16-prevfirstelem+6+1,':k')

%%
[h,p]=ttest(data(:));
[p,tbl,stats]=anova1((overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7))'.*50);
ylabel('Distance from deepest electrode (\mum)')
xlabel('Percentage current delivered')
xticklabels(fliplr({'0:100','25:75','50:50','75:25','100:0'}))
[p,tbl,stats]=anova1(overallcentroidshift_nostim'.*50);
ylabel('Distance from deepest electrode (\mum)')
xlabel('Percentage current delivered')
xticklabels(fliplr({'0:100','25:75','50:50','75:25','100:0'}))
title('Sham trials, no stimulation')
overallcentroidshift_stimonly_reshape=[overallcentroidshift_stimonly(:,1:7:counternumberofpairs*3*7),overallcentroidshift_stimonly(:,2:7:counternumberofpairs*3*7),overallcentroidshift_stimonly(:,3:7:counternumberofpairs*3*7),overallcentroidshift_stimonly(:,4:7:counternumberofpairs*3*7),overallcentroidshift_stimonly(:,5:7:counternumberofpairs*3*7),overallcentroidshift_stimonly(:,6:7:counternumberofpairs*3*7),overallcentroidshift_stimonly(:,7:7:counternumberofpairs*3*7)];


[p,tbl,stats] = anova2(overallcentroidshift_stimonly_reshape'.*50,counternumberofpairs*3);
c = multcompare(stats);

avgmove=mean(diff(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7)),'all')*-1*50;%in um
ms0=diff(overallcentroidshift(:,1:8:counternumberofpairs*3*8))*-1*50;
ms1=diff(overallcentroidshift(:,2:8:counternumberofpairs*3*8))*-1*50;
ms2=diff(overallcentroidshift(:,3:8:counternumberofpairs*3*8))*-1*50;
ms3=diff(overallcentroidshift(:,4:8:counternumberofpairs*3*8))*-1*50;
ms4=diff(overallcentroidshift(:,5:8:counternumberofpairs*3*8))*-1*50;
ms6=diff(overallcentroidshift(:,6:8:counternumberofpairs*3*8))*-1*50;
ms8=diff(overallcentroidshift(:,7:8:counternumberofpairs*3*8))*-1*50;
ms10=diff(overallcentroidshift(:,8:8:counternumberofpairs*3*8))*-1*50;
mstd0=std(ms0,[],'all')/sqrt(numel(ms0));
mstd1=std(ms1,[],'all')/sqrt(numel(ms1));
mstd2=std(ms2,[],'all')/sqrt(numel(ms2));
mstd3=std(ms3,[],'all')/sqrt(numel(ms3));
mstd4=std(ms4,[],'all')/sqrt(numel(ms4));
mstd6=std(ms6,[],'all')/sqrt(numel(ms6));
mstd8=std(ms8,[],'all')/sqrt(numel(ms8));
mstd10=std(ms10,[],'all')/sqrt(numel(ms10));

figure
hold on
bar([mean(ms0,'all'),mean(ms1,'all'),mean(ms2,'all'),mean(ms3,'all'),mean(ms4,'all'),mean(ms6,'all'),mean(ms8,'all'),mean(ms10,'all')])
er = errorbar(1:8,[mean(ms0,'all'),mean(ms1,'all'),mean(ms2,'all'),mean(ms3,'all'),mean(ms4,'all'),mean(ms6,'all'),mean(ms8,'all'),mean(ms10,'all')],[mstd0 mstd1 mstd2 mstd3 mstd4 mstd6 mstd8 mstd10],...
    [mstd0 mstd1 mstd2 mstd3 mstd4 mstd6 mstd8 mstd10]);er.Color = [0 0 0]; er.LineStyle = 'none';
ylabel('Average shift in centroid (\mum)')
xlabel('Total current delivered')
xticklabels({'','0','1','2','3','4','6','8','10'})
title('Centroid location')


figure 
hold on
bar(mean(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7).*50,2))
er = errorbar(1:5,mean(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7).*50,2),std(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7).*50,[],2)./sqrt(counternumberofpairs*3*7),...
    std(overallcentroidshift_stimonly(:,1:counternumberofpairs*3*7).*50,[],2)./sqrt(counternumberofpairs*3*7));er.Color = [0 0 0]; er.LineStyle = 'none';
ylabel('Distance from deepest electrode (\mum)')
xlabel('Percentage current delivered')
xticklabels(fliplr({'','0:100','25:75','50:50','75:25','100:0',''}))
title('Centroid location')

overallcentroidshift_stimonly=overallcentroidshift;
figure
for current=1:length(AMP)
subplot(2,4,current)
    hold on
bar(mean(overallcentroidshift_stimonly(:,current:length(AMP):counternumberofpairs*3*length(AMP)).*50,2))
er = errorbar(1:5,mean(overallcentroidshift_stimonly(:,current:length(AMP):counternumberofpairs*3*length(AMP)).*50,2),std(overallcentroidshift_stimonly(:,current:length(AMP):counternumberofpairs*3*length(AMP)).*50,[],2)./sqrt(counternumberofpairs*3),...
    std(overallcentroidshift_stimonly(:,current:length(AMP):counternumberofpairs*3*length(AMP)).*50,[],2)./sqrt(counternumberofpairs*3));er.Color = [0 0 0]; er.LineStyle = 'none';
er2(current,1:5)=std(overallcentroidshift_stimonly(:,current:length(AMP):counternumberofpairs*3*length(AMP)).*50,[],2)./sqrt(counternumberofpairs*3);
ylabel('Distance from deepest electrode (\mum)')
xlabel('Percentage current delivered')
xticks([0 1 2 3 4 5 6 7 8])
xticklabels(fliplr({'' '0:100','25:75','50:50','75:25','100:0' ''}))
title(['Centroid location: ' num2str(AMP(current)) '\muA total'])
xlim([0.5 5.5])
end
%% Across
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
overallpoolact=[];
shiftcent=[];
numrejectedchannels=0;
numacceptchan=0;
for shanknum=1:4
    for trial=1:5
        for columnnum=1:7
        check1=['S' num2str(shanknum) '_T' num2str(trial)];
        overallpoolact.(check1)=[];
        check=['T' num2str(trial)];
        savecent.(check)=[];
        end
    end
end

counternumberofpairs=0;

meananimal=[];
for ratN=6:10

        
if ratN<10
    Ratnum=['Rat_00' num2str(ratN)];
else
    Ratnum=['Rat_0' num2str(ratN)];
end
cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
D_data=dir;


for k = 3:length(D_data) % avoid using the first ones
    currD = D_data(k).name; % Get the current subdirectory name
    try
        cd([D_data(k).folder filesep currD])
        loadStimChn;
        if (stimChn(1)-stimChn(2))~=-16 || length(stimChn)~=2
            continue
        end
        load('Averagetrialresponse.mat','avgnospT')
        [normspk_all,centroidpershank]=PoolNormalisedActivity(avgnospT,startpointseconds, secondstoanalyse);
        if stimChn(1)<17
            shank=1;
        elseif stimChn(1)<33 && stimChn(1)>16
            shank=2;
            stimChn=stimChn-16;
        elseif stimChn(1)<49 && stimChn(1)>32
            shank=3;
            stimChn=stimChn-32;
        else
            shank=4;
            stimChn=stimChn-48;
        end
        
        for shanknum=1:4
            for trial=1:5
                check=['S' num2str(shanknum) '_T' num2str(trial)];
                rejectchan=find(normspk_all.(check)==0);
                numrejectedchannels=length(rejectchan)+numrejectedchannels;
                acceptchan=find(normspk_all.(check)~=0);
                numacceptchan=length(acceptchan)+numacceptchan;
                overallpoolact.(check)=[overallpoolact.(check),[zeros(16-stimChn(1),size(normspk_all.(check),2)); normspk_all.(check); zeros(32-(16-stimChn(1))-16,size(normspk_all.(check),2))]];
                meananimal.(Ratnum)(trial,shanknum)=nanmean(normspk_all.(check),'all');
            end
        end
        counternumberofpairs=counternumberofpairs+1;
        for columnnum=1:7
            for trial=1:5
                check2=['T' num2str(trial) '_I' num2str(columnnum)];
                shiftcent.(check2)=[];
            end
        end
        
        for trial=1:5
            for shanknum=1:4
                for columnnum=1:7
                    check=['S' num2str(shanknum) '_T' num2str(trial)];
                    check2=['T' num2str(trial) '_I' num2str(columnnum)];
                    shiftcent.(check2)=[shiftcent.(check2),normspk_all.(check)(:,columnnum)];
                    
                end
            end
        end
        
        for trial=1:5
            for columnnum=1:7
                check2=['T' num2str(trial) '_I' num2str(columnnum)];
                trialdat=[shiftcent.(check2)(:,1) shiftcent.(check2)(:,3) shiftcent.(check2)(:,4) shiftcent.(check2)(:,2)];
                A = trapz(trialdat,2);
                B(1:16,1)=0;
                
                for lims = 2:4
                    B(:,lims) =  trapz(trialdat(:,1:lims),2);
                end
                [~,electrodecentroid]=min(abs(B'-(A'/2)));
                check=['T' num2str(trial)];
                savecent.(check)=[savecent.(check),electrodecentroid'];
            end
        end
    catch
        
    end
    

    
    
   
end
end
shift1=savecent.T1-savecent.T2;
shift2=savecent.T2-savecent.T3;
shift3=savecent.T3-savecent.T4;
shift4=savecent.T4-savecent.T5;
meanshift=mean([shift1(shift1~=1000); shift2(shift2~=1000); shift3(shift3~=1000); shift4(shift4~=1000)]).*200;
%meanshift=[mean(savecent.T1,'all') mean(savecent.T2,'all') mean(savecent.T3,'all') mean(savecent.T4,'all') mean(savecent.T5,'all')].*200;
diffshift=diff(meanshift).*200;
stder=std(([shift1(shift1~=1000); shift2(shift2~=1000); shift3(shift3~=1000); shift4(shift4~=1000)]).*200)/sqrt(length([shift1(shift1~=1000); shift2(shift2~=1000); shift3(shift3~=1000); shift4(shift4~=1000)]));
p=ranksum(shift1(shift1~=1000),shift2(shift2~=1000));
p=ranksum(shift3(shift3~=1000),shift2(shift2~=1000));
p=ranksum(shift4(shift4~=1000),shift3(shift3~=1000));
p=ranksum(shift5(shift5~=1000),shift4(shift4~=1000));
%%
ratePool=[];
stderall=[];
SMOOTHING=2;
window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
for shanknum=1:4
    for trial=1:5
        check1=['S' num2str(shanknum) '_T' num2str(trial)];
        overallpoolact.(check1)(overallpoolact.(check1)==0)=NaN;
        overallpoolact.(check1)(15,:)=overallpoolact.(check1)(14,:);
        overallpoolact.(check1)(16,:)=overallpoolact.(check1)(14,:);
        overallpoolact.(check1)(17,:)=overallpoolact.(check1)(14,:);
         trialdata=overallpoolact.(check1);
sum((nonan_len>0 & nonan_len<51));
sum(nonan_len>0)
    nonan_len=zeros(size(trialdata,1),1);
    for i=1:size(trialdata,1)
        testlength=trialdata(i,:);
        nonan_len(i) = length(testlength(~isnan(testlength)));
        leng=50;%was 100
        nonan_len(nonan_len<leng)=0;
         if nonan_len(i)==0
             trialdata(i,:)=NaN;
            continue
         end
         input=find(~isnan(testlength));
         s = RandStream('mlfg6331_64');
        [y1,idx] = datasample(s,input,(nonan_len(i)-leng),'Replace',false);
        trialdata(i,y1)=NaN;
    end
        trPool=nanmean(overallpoolact.(check1),'all');
        lengthall = numel(overallpoolact.(check1)) - sum(isnan(overallpoolact.(check1)),'all');
        stder=nanstd(overallpoolact.(check1),0,'all')./sqrt(lengthall);
        %trPool(isnan(trPool))=[];
%         rate = conv(trPool,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
%         rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        ratePool=[ratePool,trPool];
        stderall=[stderall,stder];
    end
end

cmap=colormap('gray');
RP=[nanmean(ratePool(:,1:5));nanmean(ratePool(:,11:15));nanmean(ratePool(:,16:20));nanmean(ratePool(:,6:10))];
RP=[ratePool(1:5);ratePool(11:15);ratePool(16:20);ratePool(6:10)];
SD=[stderall(1:5);stderall(11:15);stderall(16:20);stderall(6:10)];
figure
hold on
for i=1:5
errorbar(RP(:,i),SD(:,i),'Color',cmap(i*floor((length(cmap)-50)/5),:),'marker', linS{i},'LineStyle','none')

end
for i=1:5
% Fit line to data using polyfit
c = polyfit(1:size(RP,1),RP(:,i),1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
% Evaluate fit equation using polyval
y_est = polyval(c,1:size(RP,1));
plot(y_est,'Color',cmap(i*floor((length(cmap)-50)/5),:))
end
title('Shifting centroid across shanks')
xlabel('Shank number')
ylabel('Average normalised response')
xticklabels({'1' '' '2' '' '3' '' '4'})

%%
for i=1:5
    rate1=RP(:,i);
    rate1(rate1<0)=0;
    A = trapz(1:4, rate1);
    B(1)=0;
    
    for lims = 2:4
        B(lims) =  trapz(1:lims, rate1(1:lims));
    end
    
    [~,electrodecentroid]=min(abs(B-(A/2)));
    c(i)=electrodecentroid;
end
%%

for trial=1:5
    figure
    for shanknum=1:4
        
        for ratN=6:13
            if ratN<10
                Ratnum=['Rat_00' num2str(ratN)];
            else
                Ratnum=['Rat_0' num2str(ratN)];
            end
            hold on
            try
                scatter(shanknum,meananimal.(Ratnum)(trial,shanknum));
            catch
                
            end
        end
    end
end
%%
for shanknum=1:4
    figure
    for trial=1:5
        check1=['S' num2str(shanknum) '_T' num2str(trial)];
        overallpoolact.(check1)(overallpoolact.(check1)==0)=NaN;
        overallpoolact.(check1)(15,:)=overallpoolact.(check1)(14,:);
        overallpoolact.(check1)(16,:)=overallpoolact.(check1)(14,:);
        overallpoolact.(check1)(17,:)=overallpoolact.(check1)(14,:);
         trialdata=overallpoolact.(check1);
        
        plot(trialdata')
    end
end