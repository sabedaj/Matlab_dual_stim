
%% Load in data from current directory
[electfit]=BoltzmannSigmoidFit;
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


sigtest=zeros(nChn,2);
PlotNow=0;
includeResult=1;
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
AMP_quater=[0 AMP(2:end).*(2/4)];
AMP_threequaters=[0 AMP(2:end).*(2*3/4)];

ratiosingduals=zeros(nChn,1);
ratiosingduals75=zeros(nChn,1);
ratiosingduals25=zeros(nChn,1);
ratiosingdualn=zeros(nChn,1);
ratiosingdualn75=zeros(nChn,1);
ratiosingdualn25=zeros(nChn,1);
ratiosingdualsnm=zeros(nChn,length(AMP_orig));
ratiosingduals75nm=zeros(nChn,length(AMP_orig));
ratiosingduals25nm=zeros(nChn,length(AMP_orig));
ratiosingdualnnm=zeros(nChn,length(AMP_orig));
ratiosingdualn75nm=zeros(nChn,length(AMP_orig));
ratiosingdualn25nm=zeros(nChn,length(AMP_orig));
meansig50=zeros(nChn,1);
meansig75=zeros(nChn,1);
meansig25=zeros(nChn,1);
currentavg50=zeros(nChn,length(AMP_orig));
currentavg75=zeros(nChn,length(AMP_orig));
currentavg25=zeros(nChn,length(AMP_orig));
stdersig50=zeros(nChn,1);
stdersig75=zeros(nChn,1);
stdersig25=zeros(nChn,1);
electfitratio=zeros(nChn,10);
if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
    laminar=1;
else
    laminar=0;
end
channels=find(electfit(:,1)==1);

counter=0;
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

    amp=(log(3)-electfit(elect,4))/(electfit(elect,6)+electfit(elect,5));
    %amp = ceil(amp * 2) / 2;
%     if amp*2>max(AMP)/2
%         amp=max(AMP)/2;
%     end

    if ~((amp<10)&&(amp>=0))
        amp=3;
    end
    if amp<0.5
        amp=0.5;
    end
    if amp>AMP_orig(end)/2
        amp=AMP_orig(end)/2;
    end
        [~, AMPInterestindex]=find((AMP_orig/2-amp)>=0,1); %~ is a throw away variable but it will tell you how close the value is to your chosen val

%     if AMP_orig(AMPInterestindex)>max(AMP)/2
%         AMP_orig(AMPInterestindex)=max(AMP);
%     end
%             AMP_origp=AMP_orig;
%         AMP_origp(AMPInterestindex)=-5;
    amp=(log(1/3)-electfit(elect,4))/(electfit(elect,6)+electfit(elect,5));
    
    if amp<0.5
        amp=0.5;
    end
    if ((amp<=0))
        continue
    elseif (amp>AMP_orig(end)/2)
        amp=AMP_orig(end)/2;
    end
    AMP_origp=AMP_orig;
    AMP_origp(AMPInterestindex)=-5;
        AMPInterestindex2=find((AMP_orig/2-amp)>0,1); %~ is a throw away variable but it will tell you how close the value is to your chosen val
        AMPInterestindex2=AMPInterestindex2-1;
        if isempty(AMPInterestindex2)
            AMPInterestindex2=length(AMP_orig);
        end
        %AMPInterestindex2=AMPInterestindex+1;
%     amp=AMP_orig(AMPInterestindex)*2;
%     [~, AMPInterestindexdual]=min(abs( AMP_orig-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%     AMP_origp=AMP_orig;
%         AMP_origp(AMPInterestindexdual)=-5;
%         [~, AMPInterestindex2dual]=min(abs(AMP_origp-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%     
%         if AMPInterestindexdual>AMPInterestindex2dual
%             temp=AMPInterestindex2dual;
%             AMPInterestindex2dual=AMPInterestindexdual;
%             AMPInterestindexdual=temp;
%         end
    if AMPInterestindex>AMPInterestindex2
        temp=AMPInterestindex2;
        AMPInterestindex2=AMPInterestindex;
        AMPInterestindex=temp;
    end
%     if AMPInterestindex2dual>length(electrodeampint50)
%         AMPInterestindex2dual=length(electrodeampint50);
%     end
    if AMPInterestindex2>length(electrodeampint50)
        AMPInterestindex2=length(electrodeampint50);
    end
%     AMPInterestindex2=7;
%     AMPInterestindex=6;
    % y=[electrodeampint100 electrodeampint75 electrodeampint50 electrodeampint25 electrodeampint0];
%     AMPInterestindex=2;
%     AMPInterestindex2=8;
%log((y./beta(1))./(1-(y./beta(1))))
%     stim25=log((mean(electrodeampint25(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))/(1-(mean(electrodeampint25(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))));
%     stim75=log((mean(electrodeampint75(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))/(1-(mean(electrodeampint75(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))));
%     stim50=log((mean(electrodeampint50(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))/(1-(mean(electrodeampint50(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))))+5;
%     stim25nm=log(((electrodeampint25(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))./(1-((electrodeampint25(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))));
%     stim75nm=log(((electrodeampint75(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))./(1-((electrodeampint75(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))));
%     stim50nm=log(((electrodeampint50(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))./(1-((electrodeampint50(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))))+5;
%     stim100=electrodeampint100(pos);
%     stim100nm=log(((stim100(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))./(1-((stim100(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))))+5;
%     stim100=log((mean(stim100(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))/(1-(mean(stim100(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))))+5;
%     stim0=electrodeampint0(pos);
%     stim0nm=log(((stim0(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))./(1-((stim0(AMPInterestindex:AMPInterestindex2))./electfit(elect,3))))+5;
%     stim0=log((mean(stim0(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))/(1-(mean(stim0(AMPInterestindex:AMPInterestindex2))/electfit(elect,3))))+5;
     stim25=mean(electrodeampint25(AMPInterestindex:AMPInterestindex2));

    stim75=mean(electrodeampint75(AMPInterestindex:AMPInterestindex2));
    stim50=mean(electrodeampint50(AMPInterestindex:AMPInterestindex2));
    stim25nm=(electrodeampint25(AMPInterestindex:AMPInterestindex2));
    stim75nm=(electrodeampint75(AMPInterestindex:AMPInterestindex2));
    stim50nm=(electrodeampint50(AMPInterestindex:AMPInterestindex2));
    stim100=electrodeampint100(pos);
    stim100nm=(stim100(AMPInterestindex:AMPInterestindex2));
    stim100=mean(stim100(AMPInterestindex:AMPInterestindex2));
    stim0=electrodeampint0(pos);
    stim0nm=(stim0(AMPInterestindex:AMPInterestindex2));
    stim0=mean(stim0(AMPInterestindex:AMPInterestindex2));
    
   

    
        AMP3q=AMP_orig(AMPInterestindex:AMPInterestindex2)*3/4;
        [~, Amp3qindex]=intersect(AMP,AMP3q);
        stim1003q= mean(electrodeampint100(Amp3qindex));
        stim03q=mean(electrodeampint0(Amp3qindex));
        stim1003qnm= (electrodeampint100(Amp3qindex));
        stim03qnm=(electrodeampint0(Amp3qindex));
     AMPq=AMP_orig(AMPInterestindex:AMPInterestindex2)*1/4;
    [~, Ampqindex]=intersect(AMP,AMPq);
        stim100q= mean(electrodeampint100(Ampqindex));
        stim0q=mean(electrodeampint0(Ampqindex));
        stim100qnm= (electrodeampint100(Ampqindex));
        stim0qnm=(electrodeampint0(Ampqindex));

        
        
        if AMPInterestindex~=AMPInterestindex2 && electfit(elect,1)==1 %%model
            Amprange=AMP_halfp(AMPInterestindex:AMPInterestindex2);
            Amprange_dual50=Amprange(1):0.1:Amprange(end);
            model50results(elect,1:length(Amprange_dual50))=electfit(elect,3)./(1+exp(electfit(elect,4)+electfit(elect,5).*Amprange_dual50 + electfit(elect,6).*Amprange_dual50));
            single1results(elect,1:length(Amprange_dual50))=electfit(elect,3)./(1+exp(electfit(elect,4)+electfit(elect,5).*Amprange_dual50));
            single2results(elect,1:length(Amprange_dual50))=electfit(elect,3)./(1+exp(electfit(elect,4)+electfit(elect,6).*Amprange_dual50));
            mean50modelresults(elect)=mean((model50results(elect,1:length(Amprange_dual50))-(single1results(elect,1:length(Amprange_dual50))+single2results(elect,1:length(Amprange_dual50))))./abs(single1results(elect,1:length(Amprange_dual50))+single2results(elect,1:length(Amprange_dual50))));
        elseif electfit(elect,1)==1
            Amprange_dual50=AMP_halfp(AMPInterestindex);
            model50results(elect,1:length(Amprange_dual50))=electfit(elect,3)./(1+exp(electfit(elect,4)+electfit(elect,5).*Amprange_dual50 + electfit(elect,6).*Amprange_dual50));
            single1results(elect,1:length(Amprange_dual50))=electfit(elect,3)./(1+exp(electfit(elect,4)+electfit(elect,5).*Amprange_dual50));
            single2results(elect,1:length(Amprange_dual50))=electfit(elect,3)./(1+exp(electfit(elect,4)+electfit(elect,6).*Amprange_dual50));
             mean50modelresults(elect)=mean((model50results(elect,1:length(Amprange_dual50))-(single1results(elect,1:length(Amprange_dual50))+single2results(elect,1:length(Amprange_dual50))))./abs(single1results(elect,1:length(Amprange_dual50))+single2results(elect,1:length(Amprange_dual50))));
        end
        

if (intersect(elect,channels))==elect
    counter=counter+1;
    sigtest(elect,:)=[stim50 (stim100+stim0)];
        ratiosingduals(elect)=((stim50)-(stim100+stim0))/abs(stim100+stim0);
        ratiosingduals75(elect)=(stim75-(stim1003q+stim0q))/abs(stim1003q+stim0q);
        ratiosingduals25(elect)=(stim25-(stim1003q+stim0q))/abs(stim100q+stim03q);
        
        ratiosingdualsnm(elect,1:length(stim50nm))=((stim50nm-(stim100nm+stim0nm))./abs(stim100nm+stim0nm));
        ratiosingduals75nm(elect,1:length(stim75nm))=((stim75nm-(stim1003qnm+stim0qnm))./abs(stim1003qnm+stim0qnm));
        ratiosingduals25nm(elect,1:length(stim25nm))=(stim25nm-(stim1003qnm+stim0qnm))./abs(stim100qnm+stim03qnm);
        
        meansig50(elect)=mean(((stim50nm-(stim100nm+stim0nm))./abs(stim100nm+stim0nm)));
        meansig75(elect)=mean(((stim75nm-(stim1003qnm+stim0qnm))./abs(stim1003qnm+stim0qnm)));
        meansig25(elect)=mean((stim25nm-(stim1003qnm+stim0qnm))./abs(stim100qnm+stim03qnm));
        
        stdersig50(elect)=std(((stim50nm-(stim100nm+stim0nm))./abs(stim100nm+stim0nm)))/sqrt(length(stim50nm));
        stdersig75(elect)=std(((stim75nm-(stim1003qnm+stim0qnm))./abs(stim1003qnm+stim0qnm)))/sqrt(length(stim75nm));
        stdersig25(elect)=std((stim25nm-(stim1003qnm+stim0qnm))./abs(stim100qnm+stim03qnm))/sqrt(length(stim25nm));
        
        currentavg50(elect,1:length(pos))=(electrodeampint50-(electrodeampint0(pos)+electrodeampint100(pos)))./abs(electrodeampint0(pos)+electrodeampint100(pos));
         [~,pos3]=intersect(AMP, AMP_orig.*(3/4));
         [~,pos1]=intersect(AMP, AMP_orig.*(1/4));
        currentavg75(elect,1:length(pos))=(electrodeampint75-(electrodeampint0(pos1)+electrodeampint100(pos3)))./abs(electrodeampint0(pos1)+electrodeampint100(pos3));
        currentavg25(elect,1:length(pos))=(electrodeampint25-(electrodeampint0(pos3)+electrodeampint100(pos1)))./abs(electrodeampint0(pos3)+electrodeampint100(pos1));
else
        ratiosingdualn(elect)=(stim50-(stim100+stim0))/abs(stim100+stim0);
        ratiosingdualn75(elect)=(stim75-(stim1003q+stim0q))/abs(stim1003q+stim0q);
        ratiosingdualn25(elect)=(stim25-(stim1003q+stim0q))/abs(stim100q+stim03q);
        ratiosingdualnnm(elect,1:length(stim50nm))=(stim50nm-(stim100nm+stim0nm))./abs(stim100nm+stim0nm);
        ratiosingdualn75nm(elect,1:length(stim75nm))=(stim75nm-(stim1003qnm+stim0qnm))./abs(stim1003qnm+stim0qnm);
        ratiosingdualn25nm(elect,1:length(stim25nm))=(stim25nm-(stim1003qnm+stim0qnm))./abs(stim100qnm+stim03qnm);
end

end

%in case the ratio is unreasonable (e.g. abd(stim100+stim0) is very close
%to 0)
ratiosingduals(abs(ratiosingduals)>100)=0;
ratiosingdualn(abs(ratiosingdualn)>100)=0;
ratiosingduals75(abs(ratiosingduals75)>100)=0;
ratiosingdualn75(abs(ratiosingdualn75)>100)=0;
ratiosingduals25(abs(ratiosingduals25)>100)=0;
ratiosingdualn25(abs(ratiosingdualn25)>100)=0;

ratiosingdualsnm(abs(ratiosingdualsnm)>100)=0;
ratiosingdualnnm(abs(ratiosingdualnnm)>100)=0;
ratiosingduals75nm(abs(ratiosingduals75nm)>100)=0;
ratiosingdualn75nm(abs(ratiosingdualn75nm)>100)=0;
ratiosingduals25nm(abs(ratiosingduals25nm)>100)=0;
ratiosingdualn25nm(abs(ratiosingdualn25nm)>100)=0;
save('Significant.mat','ratiosingduals','ratiosingduals75','ratiosingduals25','sigtest','counter','currentavg50','currentavg25','currentavg75','electfitratio')
save('NonSignificant.mat','ratiosingdualn','ratiosingdualn75','ratiosingdualn25')



save('Significantnm.mat','ratiosingdualsnm','ratiosingduals75nm','ratiosingduals25nm','meansig50', 'meansig75', 'meansig25','stdersig50', 'stdersig75', 'stdersig25')
save('NonSignificantnm.mat','ratiosingdualnnm','ratiosingdualn75nm','ratiosingdualn25nm')
%%
ap=ratiosingduals;
ap(isnan(ratiosingduals))=[];
norm=jbtest(ap(ap~=0));
[h,p1]=ttest(ap(ap~=0),0);



load('Significant.mat','ratiosingduals')
sigratio=ratiosingduals;
load('NonSignificant.mat','ratiosingdualn')
nonsigratio=ratiosingdualn;
nonsigratio(isnan(nonsigratio))=0;
sigratio(isnan(sigratio))=0;
nonsigratio(isinf(nonsigratio))=0;
sigratio(isinf(sigratio))=0;
electall=sum(sigratio+nonsigratio)./sum((sigratio+nonsigratio)~=0);
electsig=mean(sigratio(sigratio~=0));

electnonsig=mean(nonsigratio(nonsigratio~=0));


e1=ratiosingduals+ratiosingdualn;
for chnnn=1:nChn
    if electfit(chnnn,2)>0.75
        es75(chnnn)=e1(chnnn);
    end
    if electfit(chnnn,2)>0.8
        es8(chnnn)=e1(chnnn);
    end
    if electfit(chnnn,2)>0.85
        es85(chnnn)=e1(chnnn);
    end
    if electfit(chnnn,2)>0.9
        es9(chnnn)=e1(chnnn);
    end

end

% loadVarAmp;
% SavetoPPT=0;
% loadNREP;
% nT=n_REP_true;
% fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
% fileinfo = dir('info.rhs');
% if (datetime(fileinfo.date) < fourShank_cutoff)
%     nChn=32;
%     E_Mapnumber=0;
% else
%     E_Mapnumber=loadMapNum;
%     if E_Mapnumber>0
%         nChn=64;
%     else
%         nChn=32;
%     end
% end
% 
% TP = loadTrialParams;
% trialinfo=TrialInfo_sab;
% trialinfo(1,:)=[];
% OVERALLOOPcounter=0;
% cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
% totalnumbercond=length(trialinfo(:,2))/(2*cond); %gives the total number of condition changes i.e. 100/0 75/25 but with same electrode are classed as one condition change
% loadNORECORDELECT;
% loadStimChn;
% loadVarAmp;
% loadSingleDual;
% checkmaxsize=0;
% check=0;
% if any((NORECORDELECT+CHN(1)+1)==cell2mat(trialinfo(2:end,2)))
%     laminar=1;
% else
%     laminar=0;
% end
% 
% for count=1:length(CHN)
%     for recordcount=1:length(NORECORDELECT)
%         clear relatedtrials
%         desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)); %finds trials with desired initial electrode
%         %single stim
%         if (desiredchanneltrial(end)==length(trialinfo)) %need to remove the last trial in case there is a stim electrode pair resulting in 1 as the last position
%             desiredchanneltrial(end)=[];
%         end
%         if laminar==1
%             desiredchanneltrial_plus=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(CHN(count)+NORECORDELECT(recordcount)+1)),1);% finds triials with matching recording electrode spacing
%         else
%             desiredchanneltrial_plus=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(NORECORDELECT(recordcount))),1);% finds triials with matching recording electrode spacing
%         end
%         if (VarAmp==1)||(singanddual==1)
%             modified_desiredchanneltrial_plus= desiredchanneltrial_plus;
%             for condition=1:totalnumbercond
%                 if any(modified_desiredchanneltrial_plus<(cond*condition*2)) && (check>0)
%                     check=check+1;
%                     relatedtrials(:,check)=modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)));
%                 elseif any(modified_desiredchanneltrial_plus<(cond*condition*2))
%                     relatedtrials=modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)));
%                     check=1;
%                 end
%                 modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)))=[];
%                 
%                 if  isempty(modified_desiredchanneltrial_plus)
%                     break;
%                 end
%             end
%             if exist('relatedtrials', 'var')
%                 relatedtrials( :, all(~relatedtrials,1) ) = [];
%                 desiredchannel__singleamp=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
%                 if recordcount<2
%                     desiredchannel__singleamp(desiredchannel__singleamp<checkmaxsize)=[];
%                 end
%                 AMP=loadAMP;
%                 AMP=[-1 AMP(2:end)./2];
%                 for Amploop=1:length(AMP)%used to identify existing 75/25 amplitudes for chn 1
%                     
%                     try
%                         Desired_trial(Amploop)=(desiredchannel__singleamp(cell2mat(trialinfo(desiredchannel__singleamp,18))==AMP(Amploop))+1)/2;%array of mathcning trial number
%                     catch
%                         Desired_trial(Amploop)=0;
%                     end
%                     
%                 end
%                 desiredchannel__singleamp=Desired_trial'.*2-1;
%                 if laminar==1
%                     desiredchanneltrial=find(cell2mat(trialinfo(:,2))==(CHN(count)+NORECORDELECT(recordcount)+1)); %finds trials with desired initial electrode
%                 else
%                     desiredchanneltrial=find(cell2mat(trialinfo(:,2))==NORECORDELECT(recordcount)); %finds trials with desired initial electrode
%                 end
%                 if desiredchanneltrial(end)==length(trialinfo)
%                     desiredchanneltrial(end)=[];
%                 end
%                 desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
%                 for Amploop=1:length(AMP)%used to identify existing 75/25 amplitudes for chn 1
%                     
%                     try
%                         Desired_trial2(Amploop)=(desiredchannel__singleampmath(cell2mat(trialinfo(desiredchannel__singleampmath,18))==AMP(Amploop))+1)/2;%array of mathcning trial number
%                     catch
%                         Desired_trial2(Amploop)=0;
%                     end
%                     
%                 end
%                 desiredchannel__singleampmath=Desired_trial2'.*2-1;
% %                 if VarAmp == 1
% %                     desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<max(relatedtrials,[],'all'));
% %                 end
% %                 
% %                 if length(desiredchannel__singleampmath)~=size(relatedtrials,2)
% %                     desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
% %                     desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<(cond*2*length(relatedtrials)));
% %                 end
% 
%                 
%                 relatedtrials_all=[desiredchannel__singleampmath';relatedtrials; desiredchannel__singleamp'];
%                 relatedtrials_all=(relatedtrials_all+1)./2;
%             end
%         else
%             relatedtrials_all=(desiredchanneltrial_plus+1)./2;
%         end
% 
%         if VarAmp~=0 && ~(isempty(relatedtrials))
%             checkmaxsize=max(relatedtrials,[],'all');
%         end
%         for chn=1:nChn
%             focuschn=chn;
%             d = Depth(E_Mapnumber); chn = d(chn);
%             
%             sp = loadSpikes(chn);
%             sp = denoiseSpikes(sp,chn);
%             result_all=zeros(size(relatedtrials_all,1),size(relatedtrials_all,2));
%             %%%%%%%%%%%%%raster input
%             for counter_trials=1:size(relatedtrials_all,2)
%                 for counter_trials_row=1:size(relatedtrials_all,1)
%                     ID=relatedtrials_all(counter_trials_row,counter_trials);
%                     tID = find((cell2mat(TP(:,2)) == ID));
%                     tID(1:2:end)=[];
%                     tID=tID./2;
%                     trig = loadTrig(0);
%                     trig = trig(tID);
%                     trig(trig==-500)=[];
%                     theseTrig = trig./30;
%                     BIN = [-0 10];
%                     xdata = [];
%                     
%                     for tr = 1:length(trig)
%                         theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
%                         for i = 1:length(theseSp)
%                             xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
%                         end
%                     end
%                     if isnan(mean(xdata))
%                         result_all(counter_trials_row,counter_trials)=0;
%                     else
%                         result_all(counter_trials_row,counter_trials)=sum(xdata)/length(trig);
%                     end
%                 end
%             end
% %             if laminar==1
% %                 strsave=['Ratio' num2str(focuschn) '_stimchn' num2str(CHN(count)) '_' num2str((CHN(count)+NORECORDELECT(recordcount)+1)) '.mat'];
% %             else
% %                 strsave=['Ratio' num2str(focuschn) '_stimchn' num2str(CHN(count)) '_' num2str((NORECORDELECT(recordcount))) '.mat'];
% %             end
% %             save(strsave, 'result_all')
%             result_all=zeros(size(relatedtrials_all,1),size(relatedtrials_all,2));
%         end
%     end
% end


% 
% for chnnum=1:length(channels)
%     focuschn=channels(chnnum);
%     %find average increase around 50% point
%     AMP=loadAMP;
%     AMP=[0 AMP(2:end)];
%     amp=-electfit(focuschn,4)/(electfit(focuschn,6)+electfit(focuschn,5));
%     if (amp<10)&&(amp>=0)
%         [~, AMPInterestindex]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%         AMP(AMPInterest)=-5;
%         [~, AMPInterestindex2]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%     else
%         amp=3;
%         [~, AMPInterestindex]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%         AMP(AMPInterest)=-5;
%         [~, AMPInterestindex2]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%     end
%     if AMPInterestindex>AMPInterestindex2
%         temp=AMPInterestindex2;
%         AMPInterestindex2=AMPInterestindex;
%         AMPInterestindex=temp;
%     end
%     
%     if laminar==1
%         strsave=['Ratio' num2str(focuschn) '_stimchn' num2str(CHN(count)) '_' num2str((CHN(count)+NORECORDELECT(recordcount)+1)) '.mat'];
%     else
%         strsave=['Ratio' num2str(focuschn) '_stimchn' num2str(CHN(count)) '_' num2str((NORECORDELECT(recordcount))) '.mat'];
%     end
%     load(strsave, 'result_all')
%     stim100=mean(result_all(1,AMPInterestindex:AMPInterestindex2));
%     stim50=mean(result_all(3,AMPInterestindex:AMPInterestindex2));
%     stim0=mean(result_all(5,AMPInterestindex:AMPInterestindex2));
%     ratiosingdual(focuschn)=stim50./(stim100+stim0);
% end
% 
% save('Significant.mat','ratiosingdual')
% 
% ratiosingdual=zeros(nChn,1);
% channels=find(electfit(:,1)==0);
% for chnnum=1:length(channels)
%     focuschn=channels(chnnum);
%     %find average increase around 50% point
%     AMP=loadAMP;
%     AMP=[0 AMP(2:end)];
%     amp=-electfit(focuschn,4)/(electfit(focuschn,6)+electfit(focuschn,5));
%     if (amp<10)&&(amp>=0)
%         [~, AMPInterestindex]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%         AMP(AMPInterest)=-5;
%         [~, AMPInterestindex2]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%     else
%         amp=3;
%         [~, AMPInterestindex]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%         AMP(AMPInterest)=-5;
%         [~, AMPInterestindex2]=min(abs(AMP-amp)); %~ is a throw away variable but it will tell you how close the value is to your chosen val
%     end
%     if AMPInterestindex>AMPInterestindex2
%         temp=AMPInterestindex2;
%         AMPInterestindex2=AMPInterestindex;
%         AMPInterestindex=temp;
%     end
%     
%     if laminar==1
%         strsave=['Ratio' num2str(focuschn) '_stimchn' num2str(CHN(count)) '_' num2str((CHN(count)+NORECORDELECT(recordcount)+1)) '.mat'];
%     else
%         strsave=['Ratio' num2str(focuschn) '_stimchn' num2str(CHN(count)) '_' num2str((NORECORDELECT(recordcount))) '.mat'];
%     end
%     load(strsave, 'result_all')
%     stim100=mean(result_all(1,AMPInterestindex:AMPInterestindex2));
%     stim50=mean(result_all(3,AMPInterestindex:AMPInterestindex2));
%     stim0=mean(result_all(5,AMPInterestindex:AMPInterestindex2));
%     ratiosingdual(focuschn)=stim50./(stim100+stim0);
% end
% 
% save('NonSignificant.mat','ratiosingdual')
