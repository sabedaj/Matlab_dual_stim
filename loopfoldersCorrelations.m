folders=dir;
savecor{1}=[];
savecor{2}=[];
savecor{3}=[];
load('ElectLayerClass.mat','ElectLayerClass')
for loopfold=3:size(folders,1)
    try
        cd([folders(loopfold).folder filesep folders(loopfold).name])
    catch
        continue
    end
    
    loadAMP_all;
    trialinfo=loadTrialInfo;
    trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
    numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
    chn=unique(trialinfo(:,2));
    chn(chn==0)=[];
    count=0;
    allsingletrialIDs=[];
    IDs=cell(length(chn),1);
    for trialloop=1:numsim_elect:length(trialinfo(:,2))
        if sum(trialinfo(trialloop:trialloop+numsim_elect-1,17)~=-1)==1 || sum(trialinfo(trialloop:trialloop+numsim_elect-1,2)==0)==numsim_elect-1
            count=count+1;
            for chnloop=1:length(chn)
                if any(trialinfo(trialloop:trialloop+numsim_elect-1,2)==chn(chnloop))
                    IDs{chnloop}(count,1:3)=[trialinfo(trialloop,17) (trialloop+numsim_elect-1)./numsim_elect trialinfo(trialloop,2)];
                    if trialinfo(trialloop,17)==6
                        [maxcor,lagcor]=ratecorrelations('stim',(trialloop+numsim_elect-1)./numsim_elect);
                        sgbord=find((ElectLayerClass-1),1,'first');
                        gibord=find((ElectLayerClass-2),1,'first');
                        iwbord=find((ElectLayerClass-3),1,'first');
                        
                        savecor{ElectLayerClass(trialinfo(trialloop,2))}=savecor{ElectLayerClass(trialinfo(trialloop,2))}+maxcor;%%%%% need to align maxcor
                    end
                end
            end
            %allsingletrialIDs(count)=trialloop;%index   trialinfo(trialloop,1);
        end
    end
    
    ti=loadTrialInfo;

    

end
