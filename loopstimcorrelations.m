function stimchnall=loopstimcorrelations %lam E5 across E4
folders=dir;
%stimchnall=cell(7,1);used to find commonstimchns
for loopfold=6:size(folders,1)
    try
        cd([folders(loopfold).folder filesep folders(loopfold).name])
    catch
        continue
    end
    internalfolders=dir;
    for loopfoldinternal=3:size(internalfolders,1)
        try
            cd([internalfolders(loopfoldinternal).folder filesep internalfolders(loopfoldinternal).name])
            loadStimChn;
            stimChn_bse=stimChn;
            while any(stimChn_bse>16)
                stimChn_bse(stimChn_bse>16)=stimChn_bse(stimChn_bse>16)-16;
            end
            ti=loadTrialInfo;
            %stimchnall{loopfold-2}=[stimchnall{loopfold-2};stimChn_bse];%used to find common stimchns
        catch
            continue
        end
        if any(stimChn_bse==5) || any(stimChn_bse==4)
            chn=stimChn(stimChn_bse==5);
            if isempty(chn)
                chn=stimChn(stimChn_bse==4);
            end
            if ~isempty(chn)
                
                
                ID_CHN=cell2mat(ti(2:end,[1:2,18]));
                ids=ID_CHN((ID_CHN(:,2)==chn(1) & ID_CHN(:,3)==6),1);
                if ~isempty(ids)
                 [maxcor,lagcor]=ratecorrelations('stim',ids(1));
                 break
                end
%                if any(ID_CHN(:,2)==chn(1) & ID_CHN(:,3)==6)
%                    folders(loopfold).name
%                end
%                 index=ID_CHN(ID_CHN(:,2)==chn,1);
%                 ID_CHN(ID_CHN(:,1)==index)
                %%find ID with stimchn. then call rate correlations
            end
            
        end
        
        
    end
    chn=[];
    cd([folders(loopfold).folder filesep folders(loopfold).name])
    save('stimcorrelations.mat','maxcor')
end
end