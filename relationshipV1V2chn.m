function relationshipV1V2chn(ratestruct,savefilename,chnstoremove)
% channels in V1 are 65:128 in array while channels 1:64 are in V2 function
% is used to identify the relationship between V1 and V2
%need to identify if a particular channel in v1 always responds with the
%same V2 channel
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=cell(numfolderstotal,1);
channelsigEL=cell(numfolderstotal,1);
significantmatrixEL=cell(3,2);
v2respchn=cell(numfolderstotal,64);
E_MAP=Depth;
for numerfolders=1:numfolderstotal
    
    trialend=length(fieldnames(ratestruct{numerfolders}));
    channelsigEL{numerfolders}=nan(128,trialend);
    ratespiking{numerfolders}=nan(128,181,trialend);
    for trial=1:trialend
        trialcheck=['T' num2str(trial)];
        ratespiking{numerfolders}(:,:,trial)=nanmean(ratestruct{numerfolders}.(trialcheck)(E_MAP,:,:),3);
    end
    meancatdataEL=cat(2,mean(ratespiking{numerfolders}(:,92:136,:),2),mean(ratespiking{numerfolders}(:,137:181,:),2));
    thresh=squeeze(mean(ratespiking{numerfolders}(:,1:85,:),2)+(std(ratespiking{numerfolders}(:,1:85,:),[],2).*3));
    trialsinterest=savefilename{numerfolders}{4}(:,1);
    for chn=1:128
        
        if chn<65
            chnindex=2;
        else
            chnindex=1;
            v2respchn{numerfolders,chn-64}=zeros(65,1);
        end
        if any(chn==chnstoremove)%skip bad channels
            continue
        end
        for trialI=1:length(trialsinterest)
            trialII=trialsinterest(trialI);
            
            Datathreshold=meancatdataEL(chn,:,trialII)>thresh(chn,trialII);
            if ~isempty(Datathreshold) && sum(Datathreshold)==2%both
                channelsigEL{numerfolders}(chn,trialII)=3;
                significantmatrixEL{3,chnindex}=[significantmatrixEL{3,chnindex};ratespiking{numerfolders}(chn,:,trialII)];
            elseif ~isempty(Datathreshold) && Datathreshold(1)==1%early
                channelsigEL{numerfolders}(chn,trialII)=1;
                significantmatrixEL{1,chnindex}=[significantmatrixEL{1,chnindex};ratespiking{numerfolders}(chn,:,trialII)];
            elseif ~isempty(Datathreshold) && Datathreshold(2)==1%late
                channelsigEL{numerfolders}(chn,trialII)=2;
                significantmatrixEL{2,chnindex}=[significantmatrixEL{2,chnindex};ratespiking{numerfolders}(chn,:,trialII)];
            end
            if chn>64 && ~isnan(channelsigEL{numerfolders}(chn,trialII))
                v2respchn{numerfolders,chn-64}(1:64)=v2respchn{numerfolders,chn-64}(1:64)+~isnan(channelsigEL{numerfolders}(1:64,trialII));
                v2respchn{numerfolders,chn-64}(65,1)=v2respchn{numerfolders,chn-64}(65,1)+1;%number of times the V1 channel is responding
            end
        end
    end
    
    
    
    
    
    
    %need to look for an early-late, early-early, late-late
    %realtionship between electrodes ->the both category should be
    %included. This currently tells you whether a channel responded
    %early(1) late(2) or both(3). find how many channels in V1 respond
    %with the same channel in V2
    
    
end

%need to find the array elements in v2respchn > 1
respondingchn=cellfun(@(x) x(1:64)./x(65),v2respchn,'UniformOutput',false);
%concatenate
%plot heatmap for each folder
for numerfolders=1:numfolderstotal
    figure
    imagesc(respondingchn{numerfolders,chn})
    colormap('jet')
    colorbar
    title(['Folder ' num2str(numerfolders) ' V1-V2 relationship'])
    xlabel('V2 channel')
    ylabel('V1 channel')
    %saveas(gcf,[savefilename{numerfolders}{1} 'V1V2relationship.fig'])
    %saveas(gcf,[savefilename{numerfolders}{1} 'V1V2relationship.png'])
    %close all
end
end