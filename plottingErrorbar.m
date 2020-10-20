currentavg50(abs(currentavg50)>20)=0;
sumtoremove=sum(currentavg50,2);
currentavg50(sumtoremove==0,:)=[];


%%plotting 
load('Ratio_all.mat','meansig50array','stdersig50array')
er=mean(stdersig50array(stdersig50array~=0));
er75=mean(stdersig7525array(stdersig7525array~=0));
er=nanmean(stdersig50arrayNODUAL(stdersig50arrayNODUAL~=0));
figure

subplot(4,1,1)
hold on
title('Shank 1')
errorbar(meansig50array(1:16,:),stdersig50array(1:16,:))

subplot(4,1,4)
hold on
title('Shank 4')
errorbar(meansig50array(17:32,:),stdersig50array(17:32,:))
%legend('S1E9E15','S4E10E16','S4E5E11','S4E9E15','S4E9E15')
legend('S1E11,S4E11','S1E10,S4E10','S4E7E13','S1E9,S4E9')

subplot(4,1,2)
hold on
title('Shank 2')
errorbar(meansig50array(33:48,:),stdersig50array(33:48,:))

subplot(4,1,3)
hold on
title('Shank 3')
errorbar(meansig50array(49:64,:),stdersig50array(49:64,:))
%%
figure
%[0.628677080574876;0.119229331934484;2.21681836921736e-07;396]
%[0.667592129816793;0.0845406779344633;2.84725942482189e-14;397]
% model_series=[resultind75(1); resultind(1); resultind25(1)].*100;
% model_error = [resultind75(2); resultind(2); resultind25(2)].*100;%[0.143651201; 	0.122944593; 	0.197613506].*100;
% pval=[resultind75(3) resultind(3) resultind25(3)];%[9.54*10^(-08),	1.93*10^(-06),	0.007309218];
model_series=([resultind(1,:)']).*100;
model_error = ([resultind(2,:)']).*100;%[0.143651201; 	0.122944593; 	0.197613506].*100;
pval=[resultind(3,:)];%[9.54*10^(-08),	1.93*10^(-06),	0.007309218];
% model_series=[0.667592129816793].*100;
% model_error = [0.0845406779344633].*100;
% pval=[2.84725942482189e-14];


b=bar(model_series);

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

hold off
%%For MATLAB 2019b or later releases
hold on
% Calculate the number of bars in each group
nbars = size(model_series, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
pvall01=pval<0.01;
pvall001=pval<0.001;
pvall01(pvall001)=0;
for i=1:size(pvall01,2)
    if pvall01(:,i).*x(:,i)==0
        continue
    end
scatter(pvall01(:,i).*x(:,i), pvall01(:,i).*(model_series(i,:)+model_error(i,:)+0.05.*100),'*','k')
end
for i=1:size(pvall001,2)
    if pvall001(:,i).*x(:,i)==0
        continue
    end
scatter(pvall001(:,i).*x(:,i)-0.05, pvall001(:,i).*(model_series(i,:)+model_error(i,:)+0.001.*100),'*','k')
scatter(pvall001(:,i).*x(:,i)+0.05, pvall001(:,i).*(model_series(i,:)+model_error(i,:)+0.001.*100),'*','k')

end
hold off 
%xticklabels({'75%', '50%', '25%'})
xticklabels({'0', '1', '2', '3', '4', '5', '6', '8', '10'})
%xlabel('Ratio of current towards deep electrode')
xlabel('Current(uA)')
ylabel('Percentage increase in efficacy of dual')
%title('Ratio 25% towards favoured electrode');
ylim([-50 150])
title(['Comparison of current levels in eliciting spike rates']);
%yline(1,'--')
%title(['Shift in current effect on electric field overlap (n=', num2str(resultind75(4)) ,', ',num2str(resultind(4)),', ',num2str(resultind25(4)), ')']);
% figure 
% pval=[9.54*10^(-08),	1.93*10^(-06),	0.007309218; 3.95*10^(-07),	6.04*10^(-07),	0.005868932];
% bar(1:2,pval)
% legend('4shank Rat007','4shank Rat006','1shank Rat005')
% title('Ratio significance');
% ylabel('Pval');
% xticklabels({'Dual', 'No Dual'})


%%
figure
model_series=[0.773951234	0.729450283	0.245007022	0.673650329; 0.583327491	0.97812003	0.251271444	0.618814858];
model_error = [0.303206678	0.26077995	0.076005722	0.284282759; 0.208427275	0.312690389	0.07781031	0.259101568];

b=bar(model_series);

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end
hold off
%%For MATLAB 2019b or later releases
hold on
% Calculate the number of bars in each group
nbars = size(model_series, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');

hold off 

%xticklabels({'Dual', 'No Dual'})
xlabel('Individual trials')
title('Rat006 Ratio 50/50');

figure 
pval=[0.014077755	0.007253899	0.00227858	0.022363184; 0.007784405	0.002932655	0.002292467	0.021182887];
bar(pval)
%legend('4shank Rat007','4shank Rat006','1shank Rat005')
title('Ratio significance');
ylabel('Pval');
xlabel('Individual trials')
%xticklabels({'Dual', 'No Dual'})



%%
load('Ratio_all.mat','meansig50array','stdersig50array')
figure
errorbar(meansig50array,stdersig50array)
ylim([-2 5])
xlim([1 32])
title('Rat005')
legend('11,19','12,20','12,22','13,19','12,18','9,15')

figure
hold on
errorbar(meansig50array(:,5),stdersig50array(:,5))
errorbar(meansig50array(:,2),stdersig50array(:,2))
errorbar(meansig50array(:,3),stdersig50array(:,3))
legend('12,18','12,20','12,22')
xlim([1 32])
title('Rat005')


%%
clear resultind
meansig50array(isnan(meansig50array))=0;
meansig50array(isinf(meansig50array))=0;
meansig50array(abs(meansig50array)>20)=0;
stder=std(meansig50array(meansig50array~=0))/sqrt(length(meansig50array(meansig50array~=0)));
alllength=length(meansig50array(meansig50array~=0));
norm=jbtest(meansig50array(meansig50array~=0));
[h,p1]=ttest(meansig50array(meansig50array~=0),0);
meanaltogether=  mean(meansig50array(meansig50array~=0));
resultall=[meanaltogether; stder; p1; alllength];
for i=1:size(meansig50array,2)
x=meansig50array(:,i);
stder=std(x(x~=0))/sqrt(length(x(x~=0)));
alllength=length(x(x~=0));
norm=jbtest(x(x~=0));
[h,p1]=ttest(x(x~=0),0);
meanaltogether=  mean(x(x~=0));
resultind(:,i)=[meanaltogether; stder; p1; alllength];
end

%% currents
clear resultind
meansig50array(isnan(meansig50array))=0;
meansig50array(isinf(meansig50array))=0;
meansig50array(abs(meansig50array)>20)=0;
stder=std(meansig50array(meansig50array~=0))/sqrt(length(meansig50array(meansig50array~=0)));
alllength=length(meansig50array(meansig50array~=0));
norm=jbtest(meansig50array(meansig50array~=0));
[h,p1]=ttest(meansig50array(meansig50array~=0),0);
meanaltogether=  mean(meansig50array(meansig50array~=0));
resultall=[meanaltogether; stder; p1; alllength];
for i=1:size(meansig50array,2)
x=meansig50array(:,i);
stder=std(x(x~=0))/sqrt(length(x(x~=0)));
alllength=length(x(x~=0));
norm=jbtest(x(x~=0));
[h,p1]=ttest(x(x~=0),0);
meanaltogether=  mean(x(x~=0));
resultind(:,i)=[meanaltogether; stder; p1; alllength];
end

%%
%four shank array
Stimchnall=Stimchnall_r5;
for chnpair=1:length(Stimchnall)
    if Stimchnall(chnpair,1)<17&&Stimchnall(chnpair,2)<17
        shank=1;
    elseif Stimchnall(chnpair,1)<33&&Stimchnall(chnpair,2)<33 && Stimchnall(chnpair,1)>16&&Stimchnall(chnpair,2)>16
        shank=2;
        Stimchnall(chnpair,:)=Stimchnall(chnpair,:)-16;
    elseif Stimchnall(chnpair,1)<49&&Stimchnall(chnpair,2)<49 && Stimchnall(chnpair,1)>32&&Stimchnall(chnpair,2)>32
        shank=3;
        Stimchnall(chnpair,:)=Stimchnall(chnpair,:)-32;
    elseif Stimchnall(chnpair,1)<65&&Stimchnall(chnpair,2)<65 && Stimchnall(chnpair,1)>48&&Stimchnall(chnpair,2)>48
        shank=4;
        Stimchnall(chnpair,:)=Stimchnall(chnpair,:)-48;
    else
        %across
        Stimchnall(chnpair,:)=[0,0];
    end
end
%%

allchnpair=[normalisedstimchn_r5;normalisedstimchn_r6; normalisedstimchn_r7];
position=find((allchnpair(:,1)-allchnpair(:,2))~=-6);
allchnpair(position,1)=0;
allchnpair(position,2)=0;
shallowchannel=max(allchnpair(:,2));

chnpairnumzeros=shallowchannel-allchnpair;
chnpairnumzeros(chnpairnumzeros==shallowchannel)=-500;
r5input=r5_25;
r6input=r6_25;
r7input=r7_25;
r5_shift=zeros(128,length(Stimchnall_r5));
for i=1:length(Stimchnall_r5)
    if chnpairnumzeros(i)~=-500
r5_shift(:,i)=[zeros(chnpairnumzeros(i,2),1); r5input(:,i); zeros(128-length(r5input)-chnpairnumzeros(i,2),1)];
    end
end
r6_shift=zeros(128,length(Stimchnall_r5));
for i=1:length(Stimchnall_r6)
    if chnpairnumzeros(i+length(Stimchnall_r5))~=-500
        r6_shift(:,i)=[zeros(chnpairnumzeros(i+length(Stimchnall_r5),2),1); r6input(1:16,i);zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5),2),1);...
            zeros(chnpairnumzeros(i+length(Stimchnall_r5),2),1); r6input(17:32,i);zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5),2),1); ...
            zeros(chnpairnumzeros(i+length(Stimchnall_r5),2),1); r6input(33:48,i);zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5),2),1); ...
            zeros(chnpairnumzeros(i+length(Stimchnall_r5),2),1); r6input(49:64,i); zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5),2),1)];
    end
end
r7_shift=zeros(128,length(Stimchnall_r7));
for i=1:length(Stimchnall_r7)
    if chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6))~=-500
        r7_shift(:,i)=[zeros(chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1); r7input(1:16,i);zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1);...
            zeros(chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1); r7input(17:32,i);zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1); ...
            zeros(chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1); r7input(33:48,i);zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1); ...
            zeros(chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1); r7input(49:64,i); zeros(32-16-chnpairnumzeros(i+length(Stimchnall_r5)+length(Stimchnall_r6),2),1)];
    end
end


%
rlaminartogether=[r5_shift(1:32,:),r5_shift(33:64,:),r5_shift(65:96,:),r5_shift(97:128,:),r6_shift(1:32,:),r6_shift(33:64,:),r6_shift(65:96,:),r6_shift(97:128,:),r7_shift(1:32,:),r7_shift(33:64,:),r7_shift(65:96,:),r7_shift(97:128,:)];

sumtog=sum(rlaminartogether,1);
rlaminartogether(:,sumtog==0)=[];


%%
clear resultall
for j=1:size(rlaminartogether,1)-5
meansig50array=rlaminartogether(j:j+4,:);
meansig50array(isnan(meansig50array))=0;
meansig50array(isinf(meansig50array))=0;
meansig50array(abs(meansig50array)>100)=0;
if length(meansig50array(meansig50array~=0))<2
    resultall(:,j)=[0; 0; 0; 0];
else
stder=std(meansig50array(meansig50array~=0))/sqrt(length(meansig50array(meansig50array~=0)));
alllength=length(meansig50array(meansig50array~=0));
norm=jbtest(meansig50array(meansig50array~=0));
[h,p1]=ttest(meansig50array(meansig50array~=0),0);
meanaltogether=  mean(meansig50array(meansig50array~=0));
resultall(:,j)=[meanaltogether; stder; p1; alllength];
end

end

%%
figure
%[0.628677080574876;0.119229331934484;2.21681836921736e-07;396]
%[0.667592129816793;0.0845406779344633;2.84725942482189e-14;397]
model_series=resultall(1,:)'.*100;
model_error = resultall(2,:)'.*100;%[0.143651201; 	0.122944593; 	0.197613506].*100;
pval=[resultall(3,:)];%[9.54*10^(-08),	1.93*10^(-06),	0.007309218];
% model_series=[0.667592129816793].*100;
% model_error = [0.0845406779344633].*100;
% pval=[2.84725942482189e-14];


b=bar(model_series);

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

hold off
%%For MATLAB 2019b or later releases
hold on
% Calculate the number of bars in each group
nbars = size(model_series, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
pvall01=pval<0.05;
pvall001=pval<0.001;
pvall01(pvall001)=0;
for i=1:size(pvall01,2)
    if pvall01(:,i).*x(:,i)==0
        continue
    end
scatter(pvall01(:,i).*x(:,i), pvall01(:,i).*(model_series(i,:)+model_error(i,:)+0.05.*100),'*','k')
end
for i=1:size(pvall001,2)
    if pvall001(:,i).*x(:,i)==0
        continue
    end
scatter(pvall001(:,i).*x(:,i)-0.2, pvall001(:,i).*(model_series(i,:)+model_error(i,:)+0.05.*100),'*','k')
scatter(pvall001(:,i).*x(:,i)+0.2, pvall001(:,i).*(model_series(i,:)+model_error(i,:)+0.05.*100),'*','k')

end
hold off 
%%
hold on
scatter([14,20],[0,0], [30 30],'r', 'filled')
xlim([0.5 25.5])
xlabel('Electrode number')
ylabel('Increase of dual/single electrode (%)')
%title('Ratio 25% towards favoured electrode');
title('Electrode position and ratio comparison 50/50 trials');
% figure 
% pval=[9.54*10^(-08),	1.93*10^(-06),	0.007309218; 3.95*10^(-07),	6.04*10^(-07),	0.005868932];
% bar(1:2,pval)
% legend('4shank Rat007','4shank Rat006','1shank Rat005')
% title('Ratio significance');
% ylabel('Pval');
% xticklabels({'Dual', 'No Dual'})