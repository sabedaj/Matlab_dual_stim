
%% relate correlation variance with current steering
clear
startpointseconds=2;
secondstoanalyse=8;
folders=dir;

figure;hold on;
for sepdist=5:2:9
     checksep=['sep' num2str(sepdist)];
     varVshift.(checksep){1}=[];
     varVshift.(checksep){2}=[];
for trial=1:5
     trialcheck=['T' num2str(trial)];
    saveCsplitcurrconst.(checksep).(trialcheck)=[];
end
end
%%
folders=dir;
%%
sav_var_matrix=zeros(64,64);%[];

for loopfold=3:size(folders,1)
    try
        cd([folders(loopfold).folder filesep folders(loopfold).name])
        %load('ElectLayerClass.mat','ElectLayerClass')
        %load('VarianceCorrelation.mat','var_matrix')

    catch
        continue
    end
           load('stimcorrelations.mat','maxcor')
       
        S = spdiags(nan(64,64),0,maxcor);
 maxcor=full(S);
 var_matrix=maxcor;
    var_matrix((var_matrix==0))=nan;
    % sav_var_matrix=[sav_var_matrix; nanmean(var_matrix(1:16,:),'all') nanmean(var_matrix(17:32,:),'all') nanmean(var_matrix(33:48,:),'all') nanmean(var_matrix(49:64,:),'all')];
    sav_var_matrix(:,:,loopfold-2)=var_matrix;
%                 figure
%             surf(var_matrix)
end
sav_var_matrix=mean(sav_var_matrix,3,'omitnan');
 shankavg(:,:,1)=sav_var_matrix(1:16,1:16);
  shankavg(:,:,2)=sav_var_matrix(17:32,17:32);
 shankavg(:,:,3)=sav_var_matrix(33:48,33:48);
  shankavg(:,:,4)=sav_var_matrix(49:64,49:64);
  shankavgall=mean(shankavg,3,'omitnan');
figure; surf(shankavgall)
    %%
    folders_rat=dir;
    timesmall=3;
    FilterLength=3;
    for loopfoldrat=4:size(folders_rat,1)
        try
            cd([folders_rat(timesmall).folder filesep folders_rat(loopfoldrat).name])
            load('Averagetrialresponse.mat','avgnospT')
        catch
            continue
        end
        amp=loadAMP;
         loadStimChn;
         if ~any(amp==6) || all(abs(stimChn(1)-stimChn(2))~=[6 8 10])
             continue
         end
         sepdist=abs(stimChn(1)-stimChn(2))-1;
         checksep=['sep' num2str(sepdist)];

         [~,~,Csplit]=PoolNormalisedActivity_refactored(avgnospT,startpointseconds, secondstoanalyse, ElectLayerClass);
         peak=zeros(4,5);
         for trial=1:5
             trialcheck=['T' num2str(trial)];
             saveCsplitcurrconst.(checksep).(trialcheck)=[saveCsplitcurrconst.(checksep).(trialcheck) Csplit.C6.(trialcheck)];
              
         end
         stimChn_bse=stimChn;
            while stimChn_bse(1)>16
                stimChn_bse=stimChn_bse-16;
            end

         [~,i]=max(reshape(var_matrix(stimChn(1),:),[16,4]));
         
          [~,i2]=max(reshape(var_matrix(stimChn(2),:),[16,4]));
          diffpeaktoepos(1,1:4)=stimChn_bse(1)-i;
          diffpeaktoepos(2,1:4)=i2-stimChn_bse(2);
          %peak=nanmean(reshape((abs(diff([var_matrix(stimChn(1),:);var_matrix(stimChn(2),:)]))),[16,4]));
          varVshift.(checksep){1}=[varVshift.(checksep){1} diffpeaktoepos(1,:)+diffpeaktoepos(2,:)];%store variance in correlation chn 1    varVshift.(checksep){1}=[varVshift.(checksep){1} var_matrix(stimChn(1),:)];
              %varVshift.(checksep){2}=[varVshift.(checksep){2} var_matrix(stimChn(2),:)];%store variance in correlation chn 2
    end
%end
%%


for sepdist=5:2:9
    figure(sepdist)
    hold on
    peak=[];
    checksep=['sep' num2str(sepdist)];
for trial=1:5
     trialcheck=['T' num2str(trial)];
    smoothedenvelope=smoothData(saveCsplitcurrconst.(checksep).(trialcheck),FilterLength);
    [~,peak(:,trial)]=max(smoothedenvelope);
    
end
varVshift.(checksep){2}=(mean(diff(peak'))*50);%peak(:,3)-16-(sepdist+1)/2;%mean(diff(peak'))*50; % mean shift in peak in um per shank
scatter(varVshift.(checksep){1},varVshift.(checksep){2},'filled','k')
%          varVshift.(checksep){3}=peak(:,3)-16-(sepdist+1)/2;%mean(diff(peak'))*50; % mean shift in peak in um per shank
%          scatter(mean([varVshift.(checksep){1}; varVshift.(checksep){2}]),varVshift.(checksep){3},'filled','k')
end
%%
L=cell(32,1);
averageshift=zeros(32,1);
avg=mean([varVshift{1}; varVshift{2}]);
for i=1:32
    
    if ~isempty(varVshift{3}(avg==i/2))
        L{i}=[L{i}  varVshift{3}(avg==i/2)];
            averageshift(i)=mean(varVshift{3}(avg==i/2));
    end

end
