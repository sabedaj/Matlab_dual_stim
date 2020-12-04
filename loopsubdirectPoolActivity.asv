startpointseconds=2;
secondstoanalyse=8;

%% Laminar
overallpoolact=[];
for trial=1:5
   check1=['T' num2str(trial)];
overallpoolact.(check1)=[];
end
counternumberofpairs=0;
for ratN=6:11
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
        if (stimChn(1)-stimChn(2))~=-6
            continue
        end
        load('Averagetrialresponse.mat','avgnospT')
        [normspk_all]=PoolNormalisedActivity(avgnospT,startpointseconds, secondstoanalyse);
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
            if shanknum==shank
                continue
            end
            for trial=1:5
                check=['S' num2str(shanknum) '_T' num2str(trial)];
                check1=['T' num2str(trial)];
                overallpoolact.(check1)=[overallpoolact.(check1),[zeros(16-stimChn(1),size(normspk_all.(check),2)); normspk_all.(check); zeros(32-(16-stimChn(1))-16,size(normspk_all.(check),2))]];

            end
        end
        counternumberofpairs=counternumberofpairs+1;
    catch
        
    end
    
end

end
%%
ratePool=[];
stderall=[];
SMOOTHING=1;
window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
for trial=1:5
    check1=['T' num2str(trial)];
    overallpoolact.(check1)(overallpoolact.(check1)==0)=NaN;
    trialdata=overallpoolact.(check1);

    nonan_len=zeros(size(trialdata,1),1);
    for i=1:size(trialdata,1)
        testlength=trialdata(i,:);
        nonan_len(i) = length(testlength(~isnan(testlength)));
        nonan_len(nonan_len<100)=0;
         if nonan_len(i)==0
             trialdata(i,:)=NaN;
            continue
         end
         input=find(~isnan(testlength));
         s = RandStream('mlfg6331_64');
        [y1,idx] = datasample(s,input,(nonan_len(i)-100),'Replace',false);
        trialdata(i,y1)=NaN;
    end
    %lengthall = size(overallpoolact.(check1),2);
    

    
    trPool=nanmean(trialdata,2);
    trPool(isnan(trPool))=[];
    stder=nanstd(trialdata,0,2)./100;
    stder(isnan(stder))=[];
    rate = conv(trPool,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
    ratePool=[ratePool,rate];
    stderall=[stderall,stder];
end

cmap=colormap('hot');

figure
hold on
for i=1:5
errorbar(ratePool(:,i),stderall(:,i),'Color',cmap(i*floor((length(cmap)-90)/5),:))
end
title('Shifting centroid laminar')
xlabel('Relative electrode number')
ylabel('Average normalised response')
xline(6,':k')
xline(16,':k')
%% Across
overallpoolact=[];
for shanknum=1:4
    for trial=1:5
        check1=['S' num2str(shanknum) '_T' num2str(trial)];
        overallpoolact.(check1)=[];
    end
end
counternumberofpairs=0;
for ratN=6:11
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
        [normspk_all]=PoolNormalisedActivity(avgnospT,startpointseconds, secondstoanalyse);
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
                overallpoolact.(check)=[overallpoolact.(check),[zeros(16-stimChn(1),size(normspk_all.(check),2)); normspk_all.(check); zeros(32-(16-stimChn(1))-16,size(normspk_all.(check),2))]];

            end
        end
        counternumberofpairs=counternumberofpairs+1;
    catch
        
    end
    
end
end

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

    nonan_len=zeros(size(trialdata,1),1);
    for i=1:size(trialdata,1)
        testlength=trialdata(i,:);
        nonan_len(i) = length(testlength(~isnan(testlength)));
        nonan_len(nonan_len<100)=0;
         if nonan_len(i)==0
             trialdata(i,:)=NaN;
            continue
         end
         input=find(~isnan(testlength));
         s = RandStream('mlfg6331_64');
        [y1,idx] = datasample(s,input,(nonan_len(i)-100),'Replace',false);
        trialdata(i,y1)=NaN;
    end
        trPool=nanmean(overallpoolact.(check1),'all');
        lengthall = numel(overallpoolact.(check1)) - sum(isnan(overallpoolact.(check1)),'all');
        stder=nanstd(overallpoolact.(check1),0,'all')./lengthall;
        %trPool(isnan(trPool))=[];
%         rate = conv(trPool,window);%Used to smooth the line plots and remove volatility due to a single electrode not responding
%         rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
        ratePool=[ratePool,trPool];
        stderall=[stderall,stder];
    end
end

cmap=colormap('hot');
RP=[nanmean(ratePool(:,1:5));nanmean(ratePool(:,11:15));nanmean(ratePool(:,16:20));nanmean(ratePool(:,6:10))];
RP=[ratePool(1:5);ratePool(11:15);ratePool(16:20);ratePool(6:10)];
SD=[stderall(1:5);stderall(11:15);stderall(16:20);stderall(6:10)];
figure
hold on
for i=1:5
errorbar(RP(:,i),SD(:,i),'Color',cmap(i*floor((length(cmap)-90)/5),:))

end

title('Shifting centroid across shanks')
xlabel('Shank number')
ylabel('Average normalised response')
xticklabels({'1' '' '2' '' '3' '' '4'})