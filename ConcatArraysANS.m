%concatenate arrays
Direct_dat=dir;
rateV1all=[];
rateV2all=[];
spkarraysALLV1=[];
spkarraysALLV2=[];
for itloop=3:length(Direct_dat)
    try
    filecheck=strcmp(Direct_dat(itloop).name(end-13:end),'summarydat.mat');
    catch
        continue
    end
    if ~filecheck
        continue;
    end
    load(Direct_dat(itloop).name)
    rateV1all=[rateV1all; ratesaveV1];
    rateV2all=[rateV2all; ratesaveV2];
    
    spkarraysALLV1=[spkarraysALLV1; spkarraysaveV1];
    spkarraysALLV2=[spkarraysALLV2; spkarraysaveV2];
end
%% plotting Version 2
%firing rate
rateV1all=[rateV1all_Overlap; rateV1all_NOOverlap];
BaselineV1=mean(rateV1all(:,1:90),2);
rateV1all_bssub=rateV1all-BaselineV1;
figure
BaselineV2overlap=mean(rateV2all_Overlap(:,1:89),2);
rateV2all_bssubovelap=rateV2all_Overlap-BaselineV2overlap;
BaselineV2NOoverlap=mean(rateV2all_NOOverlap(:,1:89),2);
rateV2all_bssubNOoverlap=rateV2all_NOOverlap-BaselineV2NOoverlap;
figure
ax = axes;
hold on
ax.ColorOrder = [1 0 0; 0 1 0];

 [lineOut, ~] = stdshade(rateV2all_bssubovelap,0.2,[0 0 1],-90:90,1,ax,'across');
 [lineOut, ~] = stdshade(rateV2all_bssubNOoverlap,0.2,[0 0.75 0.75],-90:90,1,ax,'across');
line([0 0],[-5 200],'Color','k');
xlabel('Time (ms)')
ylabel('Average firing rate (Sp/s)')
legend('V2 overlapping RFs','V2 non-overlapping RFs')
xlim([-89 89])
ylim([-4 25])
set(gca,'TickDir','out');

figure
ax2 = axes;
hold on
[lineOut, ~] = stdshade(rateV1all_bssub,0.2,[1 0 0],-90:90,1,ax2,'across');
line([0 0],[-5 200],'Color','k');
xlabel('Time (ms)')
ylabel('Average firing rate (Sp/s)')
legend('V1')
xlim([-89 89])
ylim([-5 180])
set(gca,'TickDir','out');
% find response time of max if the FR passes 5SD baseline
[~,V1maxtime]=max(rateV1all_bssub,[],2);
[~,V2maxtimeOv]=max(rateV2all_bssubovelap,[],2);
[~,V2maxtimeNon]=max(rateV2all_bssubNOoverlap,[],2);
figure
ax2.ColorOrder = [1 0 0; 0 0 1];
hold on
histogram((V1maxtime-90),50,'FaceColor','r','Normalization','probability');
xlabel('Time (ms)')
ylabel('Frequency')
legend('V1')
set(gca,'TickDir','out');

figure
xlim([0 91])
ax2.ColorOrder = [1 0 0; 0 0 1];
hold on
histogram((V2maxtimeOv-90),50,'FaceColor','b','Normalization','probability');
histogram((V2maxtimeNon-90),50,'FaceColor',[0 0.75 0.75],'Normalization','probability');
xlabel('Time (ms)')
ylabel('Frequency')
legend('V2 overlapping RFs','V2 non-overlapping RFs')
set(gca,'TickDir','out');
xlim([0 91])

%%sigmoid
spkarraysALLV1=[spkarraysALLV1_NOoverlap;spkarraysALLV1_overlap];
figure
ax3=axes;
hold on
AMP=[0 2 5 6 8 10];
 [lineOut, ~] = stdshade((spkarraysALLV1-mean(spkarraysALLV1(:,1))),0.2,[1 0 0],AMP,1,ax3,'across');
[lineOut, ~] = stdshade((spkarraysALLV2_overlap-mean(spkarraysALLV2_overlap(:,1))),0.2,[0 0 1],AMP,1,ax3,'across');
 [lineOut, ~] = stdshade((spkarraysALLV2_NOoverlap-mean(spkarraysALLV2_NOoverlap(:,1))),0.2,[0 0.75 0.75],AMP,1,ax3,'across');
xlabel('Current (uA)')
ylabel('Average firing rate (Sp/s)')
legend('V1','V2 overlapping RFs','V2 non-overlapping RFs')
set(gca,'TickDir','out');
ylim([-10 75])


%% plotting
%avg raster & psth of all channels
BaselineV1=mean(rateV1all(:,1:90),2);
rateV1all_bssub=rateV1all-BaselineV1;
BaselineV2=mean(rateV2all(:,1:90),2);
rateV2all_bssub=rateV2all-BaselineV2;
figure
ax = axes;
hold on
ax.ColorOrder = [1 0 0; 0 1 0];
 [lineOut, ~] = stdshade(rateV1all_bssub,0.2,[1 0 0],-90:90,1,ax,'across');
 [lineOut, ~] = stdshade(rateV2all_bssub,0.2,[0 0 1],-90:90,1,ax,'across');
MAX=max([mean(rateV2all_bssub),mean(rateV1all_bssub)]);
line([0 0],[-5 200],'Color','k');
xlabel('Time (ms)')
ylabel('Average firing rate (Sp/s)')
legend('V1','V2')
xlim([-89 89])
ylim([-5 200])
set(gca,'TickDir','out');
% find response time of max if the FR passes 5SD baseline
[~,V1maxtime]=max(rateV1all_bssub,[],2);
[~,V2maxtime]=max(rateV2all_bssub,[],2);
figure
ax2.ColorOrder = [1 0 0; 0 0 1];
hold on
histogram((V1maxtime-90),50,'FaceColor','r','Normalization','probability');
histogram((V2maxtime-90),50,'FaceColor','b','Normalization','probability');
xlabel('Time (ms)')
ylabel('Frequency')
legend('V1','V2')
set(gca,'TickDir','out');
xlim([0 91])
%sigmoid
figure
ax3=axes;
hold on
AMP=[0 2 5 6 8 10];
 [lineOut, ~] = stdshade(spkarraysALLV1,0.2,[1 0 0],AMP,1,ax3,'across');
 [lineOut, ~] = stdshade(spkarraysALLV2,0.2,[0 0 1],AMP,1,ax3,'across');
plot(AMP,mean(spkarraysALLV1))
plot(AMP,mean(spkarraysALLV2))
xlabel('Current (uA)')
ylabel('Average firing rate (Sp/s)')
legend('V1','V2')
set(gca,'TickDir','out');
ylim([-10 100])
%Threshold

ThresholdV1=zeros(size(spkarraysALLV1,1),1);
for i=1:size(spkarraysALLV1,1)
    temp=find(fliplr(spkarraysALLV1(i,:))<=0,1,'first');
    if isempty(temp)
        ThresholdV1(i)=AMP(2);
    elseif temp~=1
        amptemp=fliplr(AMP);
        ThresholdV1(i)=amptemp(temp-1);
    else
         ThresholdV1(i)=AMP(2);
    end

end


ThresholdV2=zeros(size(spkarraysALLV2,1),1);
for i=1:size(spkarraysALLV2,1)
    temp=find(fliplr(spkarraysALLV2(i,:))<=0,1,'first');
    if isempty(temp)
        ThresholdV2(i)=AMP(2);
    elseif temp~=1
        amptemp=fliplr(AMP);
        ThresholdV2(i)=amptemp(temp-1);
    else
         ThresholdV2(i)=AMP(2);
    end

end

figure
hold on
histogram((ThresholdV1),50,'FaceColor','r','Normalization','probability');
histogram((ThresholdV2),50,'FaceColor','b','Normalization','probability');
xlabel('Current (uA)')
legend('V1','V2')
ylabel('Frequency')
text(7,0.5,['V1 N=' num2str(length(ThresholdV1)), ' V2 N=' num2str(length(ThresholdV2))])
title('Threshold current')
set(gca,'TickDir','out');

%% old threshold with 50% sigmoid
halfmaxV1=max(spkarraysALLV1,[],2)/2;
[~,V1_I]=min(abs(spkarraysALLV1-halfmaxV1),[],2);
AMP=[0 2 5 6 8 10];
thresholdsV1=AMP(V1_I);
CountThresholdsV1=zeros(1,length(AMP));
for i=1:length(AMP)
CountThresholdsV1(i)=sum(thresholdsV1==AMP(i));
end

figure
bar(AMP,CountThresholdsV1)
xlabel('Current (uA)')
halfmaxV2=max(spkarraysALLV2,[],2)/2;
[~,V2_I]=min(abs(spkarraysALLV2-halfmaxV2),[],2);
thresholdsV2=AMP(V2_I);
CountThresholdsV2=zeros(1,length(AMP));
for i=1:length(AMP)
CountThresholdsV2(i)=sum(thresholdsV2==AMP(i));
end
hold on
bar(AMP,CountThresholdsV2./sum(CountThresholdsV2))
xlabel('Current (uA)')
set(gca,'TickDir','out');