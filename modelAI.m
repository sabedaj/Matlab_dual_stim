

[electfit]=BoltzmannFit_v2; % base off single electrode stimualtion due to non-symmetrical responses
%%
electrodeampint100all=[];
electrodeampint75all=[];
electrodeampint50all=[];
electrodeampint25all=[];
electrodeampint0all=[];
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
            nChn=64;
            loadCHN;
            chosenstimchn=1;
            loadNORECORDELECT;
            load('truedatastruct.mat','truedatastruct')
            step=1;
            [electfit]=BoltzmannFit_v2; % base off single electrode stimualtion due to non-symmetrical responses
            for elect=1:nChn
                if electfit(elect,1)==1
                    check=['T100_' num2str(CHN(chosenstimchn))];
                    electrodeampint100all=[electrodeampint100all; truedatastruct.(check)(elect,1:step:end)];
                    check=['T75_25_' num2str(CHN(chosenstimchn))];
                    electrodeampint75all=[electrodeampint75all; truedatastruct.(check)(elect,1:step:end)];
                    check=['T25_75_'  num2str(CHN(chosenstimchn))];
                    electrodeampint25all=[electrodeampint25all,; truedatastruct.(check)(elect,1:step:end)];
                    check=['T50_50_'  num2str(CHN(chosenstimchn))];
                    electrodeampint50all=[electrodeampint50all; truedatastruct.(check)(elect,1:step:end)];
                    check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
                    electrodeampint0all=[electrodeampint0all; truedatastruct.(check)(elect,1:step:end)];
                end
            end
        catch
        end
    end
end
%%
load('truedatastruct.mat','truedatastruct')
filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
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
AMP=loadAMP;
loadAMP_all;
loadCHN;
chosenstimchn=1;
trialinfo=loadTrialInfo;
AMP_orig=[0 AMP(2:end)];
AMP=[0 AMP_all(2:end)'];
AMP_double=[0 AMP(2:end).*2];
AMP_half=AMP_double./2;


currentavg50=zeros(nChn,length(AMP_orig));
currentavg75=zeros(nChn,length(AMP_orig));
currentavg25=zeros(nChn,length(AMP_orig));

if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
    laminar=1;
else
    laminar=0;
end



%%
channels=find(electfit(:,1)==1);
for elect=1:nChn
    check=['T100_' num2str(CHN(chosenstimchn))];
    electrodeampint100=truedatastruct.(check)(elect,1:end);
    check=['T75_25_' num2str(CHN(chosenstimchn))];
    electrodeampint75=truedatastruct.(check)(elect,1:end);
    check=['T25_75_'  num2str(CHN(chosenstimchn))];
    electrodeampint25=truedatastruct.(check)(elect,1:end);
    check=['T50_50_'  num2str(CHN(chosenstimchn))];
    electrodeampint50=truedatastruct.(check)(elect,1:end);
    if laminar==1
        check=['T100_' num2str(CHN(chosenstimchn)+NORECORDELECT(1)+1)];
    else
        check=['T100_' num2str(NORECORDELECT(1))];
    end
    electrodeampint0=truedatastruct.(check)(elect,1:end);
    [~,pos]=intersect(AMP, AMP_orig./2);
    AMP_halfp=AMP_half(pos);

        %if (intersect(elect,channels))==elect
            currentavg50(elect,1:length(pos))=(electrodeampint50-(electrodeampint0(pos)+electrodeampint100(pos)))./abs(electrodeampint0(pos)+electrodeampint100(pos));
            [~,pos3]=intersect(AMP, AMP_orig.*(3/4));
            [~,pos1]=intersect(AMP, AMP_orig.*(1/4));
            currentavg75(elect,1:length(pos))=(electrodeampint75-(electrodeampint0(pos1)+electrodeampint100(pos3)))./abs(electrodeampint0(pos1)+electrodeampint100(pos3));
            currentavg25(elect,1:length(pos))=(electrodeampint25-(electrodeampint0(pos3)+electrodeampint100(pos1)))./abs(electrodeampint0(pos3)+electrodeampint100(pos1));
        %end
        
    eFit=elect;
    %if electfit(eFit,1)==1
%         x1fit = 0:0.00005:0.01;%voltage in volts
%         x2fit = 0:0.00005:0.01;
        
        x1fit = 0:0.1:2;%voltage in volts
        x2fit = 0:0.1:2;
        [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
        beta=electfit(eFit,3:8);
       %YFIT = beta(1)./(1+exp(beta(2)+beta(3).*X1FIT + beta(4).*X2FIT));
         YFIT =beta(1)./(1+exp(beta(2)+(beta(3).*X1FIT)+(beta(4).*X2FIT)+beta(5).*((X1FIT)).*(X2FIT)));%beta(6).*((((10-X1FIT).^2).*((10-X2FIT).^2)));%
        %YFIT = beta(1)./(1+beta(4).*exp(beta(2).*X1FIT+beta(3).*X2FIT)).^2;
       % YFIT = beta(1).*((X1FIT.^beta(2))+(X2FIT.^beta(3)))./(beta(4).^(beta(5).*beta(6))+((X1FIT.^(beta(2).*beta(7)))+(X2FIT.^(beta(3).*beta(8)))));

        sing1 = YFIT(1,:);
        sing2 = YFIT(:,1);
        dual=zeros(1,length(YFIT));
        for i=1:length(YFIT)
            dual(i)=YFIT(i,i);
        end
        
        AIpost= (dual - (sing1 + sing2)) ./ (sing1 + sing2);
         AIpost= (YFIT - (sing1 + sing2)) ./ abs(sing1 + sing2);
        figure
        mesh(X1FIT,X2FIT,AIpost);
        hold on
        scatter3(AMP_orig(2:end)./2,AMP_orig(2:end)./2,currentavg50(eFit,2:end),'filled')
        scatter3(AMP_orig(2:end)*3/4,AMP_orig(2:end)*1/4,currentavg75(eFit,2:end),'filled')
        scatter3(AMP_orig(2:end)*1/4,AMP_orig(2:end)*3/4,currentavg25(eFit,2:end),'filled')
        title(num2str(elect))
                figure
        mesh(X1FIT,X2FIT,(YFIT));
        hold on
        scatter3(AMP_halfp(2:end),AMP_halfp(2:end),electrodeampint50(2:end),'filled')
        scatter3(AMP(2:end),zeros(1,length(AMP(2:end))),electrodeampint100(2:end),'filled')
        scatter3(zeros(1,length(AMP(2:end))),AMP(2:end),electrodeampint0(2:end),'filled')
        scatter3(AMP(pos3(2:end)),AMP(pos1(2:end)),electrodeampint75(2:end),'filled')
        scatter3(AMP(pos1(2:end)),AMP(pos3(2:end)),electrodeampint25(2:end),'filled')
        title(num2str(elect))
    %end
    
end
