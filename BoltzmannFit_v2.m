function [electfit]=BoltzmannFit_v2
warning('on','MATLAB:nearlySingularMatrix')
warning('on','stats:nlinfit:IllConditionedJacobian')
warning('on','stats:nlinfit:IterationLimitExceeded')

PlotNow=0;
includeResult=1;
resultonly=0;
trialinfo=loadTrialInfo;
trialinfo(1,:)=[];
loadCHN;
loadStimChn;
chosenstimchn=1;
filepath = pwd;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir([filepath, filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=loadMapNum;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
loadNORECORDELECT;
if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
    laminar=1;
else
    laminar=0;
end


%perform model fit

AMP=loadAMP;
loadAMP_all;
loadCHN;
trialinfo=loadTrialInfo;
AMP_orig=[0 AMP(2:end)];
AMP=[0 AMP_all(2:end)'];
AMP_double=[0 AMP(2:end).*2];
AMP_half=AMP_double./2;
AMP_quater=[0 AMP(2:end).*(2/4)];
AMP_threequaters=[0 AMP(2:end).*(2*3/4)];
load('truedatastruct.mat','truedatastruct')
MSEtotal=0;
electbeta=zeros(nChn,5);
electfit=zeros(nChn,8);
step=1;
[~,pos]=intersect(AMP_double, AMP_orig);
AMP_threequatersp=AMP_threequaters(pos);
AMP_quaterp=AMP_quater(pos);
AMP_halfp=AMP_half(pos);


if includeResult==1
    x1=[AMP(2:step:end) AMP_threequatersp(2:step:end) AMP_halfp(2:step:end) AMP_quaterp(2:step:end) zeros(1,length(AMP(2:step:end)))];
    x2=[zeros(1,length(AMP(2:step:end))) AMP_quaterp(2:step:end) AMP_halfp(2:step:end) AMP_threequatersp(2:step:end) AMP(2:step:end)];
    X=[x1' x2'];
elseif resultonly==1
    x1=[AMP_threequatersp(2:step:end) AMP_halfp(2:step:end) AMP_quaterp(2:step:end)];
    x2=[AMP_quaterp(2:step:end) AMP_halfp(2:step:end) AMP_threequatersp(2:step:end)];
    X=[x1' x2'];
else
    x1=[AMP(2:end) zeros(1,length(AMP(2:end)))];
    x2=[zeros(1,length(AMP(2:end))) AMP(2:end)];
    X=[x1' x2'];
end
loadStimChn;
stimchns=stimChn';

for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17&&stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
            shank=4;
            stimchns(i,:)=stimchns(i,:)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
            shank=2;
            stimchns(i,:)=stimchns(i,:)-32;
        elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
            shank=3;
            stimchns(i,:)=stimchns(i,:)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceEs=[distanceEs(:,1),distanceEs(:,4),distanceEs(:,2),distanceEs(:,3)];%change to electrode mapping
        distanceEd=[distanceEd(:,1),distanceEd(:,4),distanceEd(:,2),distanceEd(:,3)];%change to electrode mapping
        for currentsi=1:length(x1)
            currents1=x1(currentsi)*10^-6;
            currents2=x2(currentsi)*10^-6;
            check=['Es' num2str(i) '_' num2str(currents2/(10^-6))];
            check=strrep(check,'.','_');
            distpair.(check)=((currents2)./(distanceEs*10^-6*0.3*2*pi));
            check=['Ed' num2str(i) '_' num2str(currents1/(10^-6))];
            check=strrep(check,'.','_');
            distpair.(check)=((currents1)./(distanceEd*10^-6*0.3*2*pi));
        end
    end
end

for elect=1:nChn
    check=['T100_' num2str(CHN(chosenstimchn))];
    electrodeampint100=truedatastruct.(check)(elect,2:step:end);
    check=['T75_25_' num2str(CHN(chosenstimchn))];
    electrodeampint75=truedatastruct.(check)(elect,2:step:end);
    check=['T25_75_'  num2str(CHN(chosenstimchn))];
    electrodeampint25=truedatastruct.(check)(elect,2:step:end);
    check=['T50_50_'  num2str(CHN(chosenstimchn))];
    electrodeampint50=truedatastruct.(check)(elect,2:step:end);
    if laminar==1
        check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
    else
        check=['T100_' num2str(NORECORDELECT(1))];
    end
    electrodeampint0=truedatastruct.(check)(elect,2:step:end);
    
    if includeResult==1
        y=[electrodeampint100 electrodeampint75 electrodeampint50 electrodeampint25 electrodeampint0];
    elseif resultonly==1
        y=[electrodeampint75 electrodeampint50 electrodeampint25];
    else
        y=[electrodeampint100 electrodeampint0];
    end
    
    V1=zeros(1,length(x2));
    V2=zeros(1,length(x1));
    for currents=1:length(x1)
        check=['Es' num2str(i) '_' num2str(x2(currents))];
        check=strrep(check,'.','_');
        V1(currents)=distpair.(check)(elect);
        check=['Ed' num2str(i) '_' num2str(x1(currents))];
        check=strrep(check,'.','_');
        V2(currents)=distpair.(check)(elect);
    end
    X=[V2' V1'];
    X=[x1' x2'];
     %modelfun = @(b,X) b(1)./(1+exp(b(2)+(b(3).*X(:,1)+b(4).*X(:,2)))); %SIGMOIDAL  b1 is asymptote, b2=x50/slope b3=-1/slope b2=-1/slope  This is boltzmann sigmoid
    modelfun = @(b,X) b(1)./(1+exp(b(2)+(b(3).*X(:,1)+b(4).*X(:,2)+b(5).*(X(:,1)).*(X(:,2)))));%+b(5).*(10-X(:,1)).*(10-X(:,2));%+b(5).*((((10-X(:,1)).^2).*((10-X(:,2)).^2))); %SIGMOIDAL  b1 is asymptote, b2=x50/slope b3=-1/slope b2=-1/slope  This is boltzmann sigmoid
    %modelfun = @(b,X) b(1)-b(1)*(exp(b(2).*X(:,1))+exp(b(3).*X(:,2)));%+b(4).*(X(:,1).*X(:,2)))); 
    %modelfun = @(b,X) b(1)./(1+b(4).*exp(b(2).*X(:,1)+b(3).*X(:,2))).^2;
    %modelfun = @(b,X) b(1)+b(2).*X(:,1)+b(3).*X(:,2)+b(4).*X(:,1).^2+b(5).*X(:,2).^2;%quadratic

     %modelfun = @(b,X) b(1).*((X(:,1).^b(2))+(X(:,2).^b(3)))./(b(4).^(b(5).*b(6))+((X(:,1).^(b(2).*b(7)))+(X(:,2).^(b(3).*b(8)))));
    %%Naka rushton modified^  can remove 6 7 8 parameters for original that
    %%models the data
%modelfun = @(b,X) b(1).*((X(:,1)+X(:,2)). ^b(2))./(b(3).^(b(2))+((X(:,1)+X(:,2)).^b(2))); %Naka rushton original

     %beta0 = [350 2 2 0.008 2 1 1 1]; %Naka rushton modified
    %beta0 = [350 2 0.008]; %Naka original
    %beta0 = [350 2 2 4 2 1 1 1];
    beta0 = [350 1 1 1 1];
    try
        mdl = fitnlm(X,y',modelfun,beta0);
    catch
        continue
    end
    fprintf([num2str(elect) '\n']);
    beta=mdl.Coefficients.Estimate;
    electbeta(elect,:)=beta;
    pval = coefTest(mdl);
    rsq = mdl.Rsquared.Adjusted;
    
    if PlotNow==1
        x1fit = min(V2):0.00002:max(V2);
        x2fit = min(V1):0.00002:max(V1);
        x1fit = 0:0.25:10;
        x2fit = 0:0.25:10;
        
        [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
        %YFIT =  beta(1)+beta(2).*X1FIT+beta(3).*X2FIT+(beta(4).*X1FIT.^2) + (beta(5).*X2FIT.^2); %quadratic

        %YFIT =  beta(1).*((X1FIT.^beta(2))+(X2FIT.^beta(3)))./(beta(4).^(beta(5))+(X1FIT.^(beta(2)*beta(7)))+(X2FIT.^(beta(3)*beta(8))));
        %YFIT = beta(1).*((X1FIT+X2FIT).^beta(2))./(beta(3).^(beta(2))+((X1FIT+X2FIT).^beta(2)));
        %YFIT = beta(1)./(1+exp(beta(2)+beta(3).*X1FIT + beta(4).*X2FIT));
        %YFIT =beta(1)./(1+exp(beta(2)+(beta(3).*X1FIT)+(beta(4).*X2FIT)+(beta(5).*(X1FIT.*X2FIT))));
        %YFIT = beta(1)-beta(1)*exp((beta(3).*X1FIT+beta(4).*X2FIT)+beta(5).*(X1FIT.*X2FIT)); 
        %YFIT = beta(1)./(1+beta(4).*exp(beta(2).*X1FIT+beta(3).*X2FIT)).^2;
        YFIT = beta(1).*((X1FIT.^beta(2))+(X2FIT.^beta(3)))./(beta(4).^(beta(5).*beta(6))+((X1FIT.^(beta(2).*beta(7)))+(X2FIT.^(beta(3).*beta(8)))));

        figure
        mesh(X1FIT,X2FIT,(YFIT));
        hold on
        scatter3(x1,x2,y,'filled','b')
        scatter3(V2,V1,y,'filled','b')
        %scatter3(x1,x2,log((y./beta(1))./(1-(y./beta(1)))),'filled')
        %xlabel(['Channel ' num2str(CHN(chosenstimchn)) ' Current level (uA)'])
        ylabel(['Channel ' num2str(CHN(chosenstimchn)) '  Voltage (V)'])
        %ylabel(['Channel ' num2str(stimChn(chosenstimchn+1)) ' Current level (uA)'])
        xlabel(['Channel ' num2str(stimChn(chosenstimchn+1)) ' Voltage (V)'])
        zlabel('Sp/s')
        title(['Spike rates from electrode ' num2str(elect)])
        
        
        figure
        hold on
        x1fit = AMP_orig;
        x2fit = AMP_orig;
        YFIT = beta(1)./(1+exp(beta(2)+beta(3).*x1fit + beta(4).*x2fit));
        YFIT1 = beta(1)./(1+exp(beta(2)+beta(3).*x1fit));
        YFIT2 = beta(1)./(1+exp(beta(2) + beta(4).*x2fit));
        plot(AMP_orig,YFIT, 'r')
        plot(AMP_orig,YFIT1,'k-.')
        plot(AMP_orig,YFIT2,'k--')
        plot(AMP_orig,YFIT2+YFIT1,'k')
        xlabel('current (uA)')
        ylabel('Sp/s')
        title(['Model spike rates from electrode ' num2str(elect)])
        plot(AMP_halfp(1:step:end) ,electrodeampint50, 'r:')
        plot(AMP_halfp(1:step:end),electrodeampint0(pos)+electrodeampint100(pos),'k:')
        
        scatter(AMP_halfp(1:step:end) ,electrodeampint50,10, 'r')
        scatter(AMP_halfp(1:step:end),electrodeampint0(pos)+electrodeampint100(pos),10,'k')
        
        legend('Dual model', ['E' num2str(stimChn(chosenstimchn)) ' model'],['E' num2str(stimChn(chosenstimchn+1)) 'model'],['E' num2str(stimChn(chosenstimchn)) '+E' num2str(stimChn(chosenstimchn+1)) ' model'],'Dual real',['E' num2str(stimChn(chosenstimchn)) '+E' num2str(stimChn(chosenstimchn+1)) ' real'])
        
        
    end
    
     [warnmsg, msgid] = lastwarn;
    
    if  strcmp(msgid,'')  && (rsq>0.5) && (beta(1)<10000) && ((-beta(2)/beta(3))>0) && ((-beta(2)/beta(4))>0)  && ((beta(3)/beta(4))>(1/5) && (beta(3)/beta(4))<5)% (rsq>0.7)  &&&& (pval<0.05) %&& ((max(electrodeampint0)/max(electrodeampint100))<5 && (max(electrodeampint0)/max(electrodeampint100))>0.2) %&& ((beta(3)/beta(4))>0.5 && (beta(3)/beta(4))<2) % ((max(electrodeampint0)/max(electrodeampint100))<3 && (max(electrodeampint0)/max(electrodeampint100))>0.333)
        electfit(elect,1)=1;
        fprintf('Good fit \n')
        
    end
    lastwarn(''); % reset warning state
    electfit(elect,2)=rsq;
    electfit(elect,3:7)=beta;
end

save('Boltzsig.mat','electfit')
end

