sp=loadSpikes;
order=Depth(1);
spsort=sp(order);
%%
trialnum=1:25;
chnavg=zeros(64,64);
%%
ID=10;
trig = loadTrig(0);
TP = loadTrialParams;
tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
theseTrig = trig(tID)./30;
nT=length(theseTrig);


for chnoverall=1:64
spchn1=spsort{chnoverall};
spchn1(:,1)=round(spchn1(:,1));
val=zeros(64,20000/100);
pos=zeros(64,20000/100);
for lag=1:nT
    spind1=zeros(100,1);
    if lag~=1
        spind1(ceil(round(spchn1(spchn1(:,1)>(theseTrig(lag)) & spchn1(:,1)<=theseTrig(lag)-1+100,1),0)-(theseTrig(lag))))=1;
    else
        spind1(ceil(round(spchn1(spchn1(:,1)>(theseTrig(lag)) & spchn1(:,1)<(theseTrig(lag)-1)+100,1),0)-(theseTrig(lag))))=1;
    end
    for chns=1:64
        spind2=zeros(100,1);
        spchn2=spsort{chns};
        spchn2(:,1)=round(spchn2(:,1));
        if lag~=1
            spind2(ceil(round(spchn2(spchn2(:,1)>(theseTrig(lag)) & spchn2(:,1)<=theseTrig(lag)-1+100,1),0)-(theseTrig(lag))))=1;
        else
            spind2(ceil(round(spchn2(spchn2(:,1)>(theseTrig(lag)) & spchn2(:,1)<(theseTrig(lag)-1)+100,1),0)-(theseTrig(lag))))=1;
        end
        if sum(spind1)~=0 && sum(spind2)~=0
            [xc,lags] = xcorr(spind1,spind2,'coeff');
            %stem(lags,xc)
            val(chns,lag)=max(xc);
            if ~isnan(max(xc))
                pos(chns,lag)=find(xc==max(xc),1,'first');
            end
        end
    end
%     if sum(sum(val,1)~=0)>=trialnum(end)
%         break
%     end
end
pos(:,sum(val,1)==0)=[];
val(:,sum(val,1)==0)=[];

chnavg(chnoverall,:)=mean(val,2);
end
%% artificial matrix
chnavg_artificial=zeros(64,64);
for i=1:16%17:32
chnavg_artificial(i,1:16)=1;
end
for i=17:32
chnavg_artificial(i,17:32)=1;
end
for i=33:48
chnavg_artificial(i,33:48)=1;
end
for i=49:64
chnavg_artificial(i,49:64)=1;
end
figure;
surf(1:64,1:64,chnavg_artificial)
ylim([1 64])
xlim([1 64])
view(0,90)
xlabel('channel #')
ylabel('channel #')
axis square
colorbar
%%
distchn=zeros(64,64);
for focuschn=1:64
    if focuschn>16 && focuschn<33
        shank=4;
        chn=focuschn-16;
    elseif focuschn>32 && focuschn<49
        shank=2;
        chn=focuschn-32;
    elseif focuschn>48 && focuschn<65
        shank=3;
        chn=focuschn-48;
    else
        shank=1;
        chn=focuschn;
    end
    for secondchn=1:64
        if secondchn>16 && secondchn<33
            shank2=4;
            chn2=secondchn-16;
        elseif secondchn>32 && secondchn<49
            shank2=2;
            chn2=secondchn-32;
        elseif secondchn>48 && secondchn<65
            shank2=3;
            chn2=secondchn-48;
        else
            shank2=1;
            chn2=secondchn;
        end
        distchn(focuschn,secondchn)=sqrt((abs(chn-chn2)*50)^2+(abs(shank-shank2)*200)^2);
    end
end
normdistchn=distchn./max(distchn,[],'all');
figure; surf(1:64,1:64,normdistchn);
[A,X] = sort(distchn(:));
B=chnavg(:);
B=B(X);
figure %distance vector
scatter(A,B)
ylabel('Correlation')
xlabel('Distance (\mum)')

%%
chnavg_artificial=zeros(64,64);
for i=1:64
    chnavg_artificial=0;
end
%%
%load("ElectLayerClass.mat")
%ElectLayerClass(chnavg)
figure; surf(1:64,1:16,chnavg(1:16,1:16))
ylim([1 16])
xlim([1 16])
view(0,90)
xlabel('channel #')
ylabel('channel #')
yline(6,'r',LineWidth=2)
xline(6,'r',LineWidth=2)
yline(10,'r',LineWidth=2)
xline(10,'r',LineWidth=2)
axis square
colorbar

%%
figure; 
trialnum=1:100; %size(val,2)
subplot(1,4,1)
surf(trialnum,1:16,val(1:16,trialnum))  
ylim([1 16])
 view(0,90)
 ylabel('Electrode #')
xlabel('Trial #')
caxis([0 1])
subplot(1,4,4)
surf(trialnum,1:16,val(17:32,trialnum))  
 view(0,90)
  ylim([1 16])
xlabel('Trial #')
caxis([0 1])
subplot(1,4,2)
surf(trialnum,1:16,val(33:48,trialnum))   
 view(0,90)
  ylim([1 16])
xlabel('Trial #')
caxis([0 1])
subplot(1,4,3)
surf(trialnum,1:16,val(49:64,trialnum))
 view(0,90)
 ylim([1 16])
xlabel('Trial #')
caxis([0 1])

%%
figure
trialnum=1:100; %size(val,2)
subplot(1,4,1)
surf(trialnum,1:16,pos(1:16,trialnum))  
ylim([1 16])
 view(0,90)
 ylabel('Electrode #')
xlabel('Trial #')
caxis([0 200])
subplot(1,4,4)
surf(trialnum,1:16,pos(17:32,trialnum))  
 view(0,90)
  ylim([1 16])
xlabel('Trial #')
caxis([0 200])
subplot(1,4,2)
surf(trialnum,1:16,pos(33:48,trialnum))   
 view(0,90)
  ylim([1 16])
xlabel('Trial #')
caxis([0 200])
subplot(1,4,3)
surf(trialnum,1:16,pos(49:64,trialnum))
 view(0,90)
 ylim([1 16])
xlabel('Trial #')
caxis([0 200])


%%
%%
for i=1:64
    chnavg(i,:)=chnavg(:,i);
end