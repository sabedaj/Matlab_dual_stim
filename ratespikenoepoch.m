function ratespikenoepoch(ratestruct,savefilename,chnrange,excitesupress, chnexclude,normalisedat)
% used to plot the average spike rate of a channel regardless of epoch
% input: ratestruct, savefilename, chnrange, excitesupress,
% chnexclude,normalisedat output: none
if normalisedat==1
    multiplyspk=1;
else
    multiplyspk=1000;
end
D_data=dir;
if length(D_data)>size(ratestruct)
    cd([D_data(12).folder filesep D_data(12).name;])
end
E_MAP=Depth(1);
alldata=cell(length(savefilename{15}{2}.AMP),1);
groupdata=cell(length(savefilename{15}{2}.AMP),1);
correlationdata=cell(length(savefilename{15}{2}.AMP),1);
stimchnsignificant=cell(length(savefilename{15}{2}.AMP),1);
allV1electsignificantstimchn=cell(length(savefilename{15}{2}.AMP),1);
alldata10uA=cell(5,1);
totalchncount=zeros(5,1);
heatmap_centroid=cell(5,1);
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=cell(numfolderstotal,1);
for ampit=1:length(savefilename{15}{2}.AMP)
    heatmap_centroid{ampit}=nan(31,7,500);
    ampinterest=savefilename{15}{2}.AMP(ampit);
    groupdata{ampit}=cell(9,1);
    correlationdata{ampit}=cell(9,1);
    stimchnsignificant{ampit}=cell(9,1);
    allV1electsignificantstimchn{ampit}=cell(9,1);
    alldata10uA{ampit}=cell(9,1);
    
    for numerfolders=1:numfolderstotal
        
        trialend=length(fieldnames(ratestruct{numerfolders}));
        for trial=1:trialend
            trialcheck=['T' num2str(trial)];
            ratespiking{numerfolders}(:,:,trial)=nanmean(ratestruct{numerfolders}.(trialcheck),3);
        end
        [m,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest)); ampinterest=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp
        if excitesupress==1
            %determines if spiking increases with increases in current overall
            %- just has to have a net increase between 2 and 10
            rateecurrent=nan(128,length(savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1)),5);
            for ampit1=1:length(savefilename{15}{2}.AMP)
                ampinterest1=savefilename{15}{2}.AMP(ampit1);
                [m,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest1)); ampinterest1=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp
                
                rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest1),1));
                rateecurrent(1:size(rateAMPua,1),:,ampit1)=mean(rateAMPua(:,92:185,:),2);
            end
            validchannels=sum(diff(rateecurrent,[],3),3,'omitnan')>0;
        end
        rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1));
        rateAMPua(chnexclude,:,:)=nan;
        if normalisedat==1
            rateAMPua=rateAMPua./max(squeeze(max(rateAMPua,[],2)),[],2);
        end
        simchnall=savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),3);
        stimchnall_notsorted=E_MAP(simchnall);
        totalchncount(ampit)=length(chnrange)*size(rateAMPua,3)+totalchncount(ampit);
        meancatdata=mean(rateAMPua(:,92:185,:),2);
        sig=false(128,size(meancatdata,3));
        
        if excitesupress==0
            thresh=squeeze(mean(rateAMPua(:,1:85,:),2)-(std(rateAMPua(:,1:85,:),[],2).*1)); %%%%%NOTE LOWER THRESHOLD
            sig=(meancatdata)<thresh;
        else
            thresh=squeeze(mean(rateAMPua(:,1:85,:),2)+(std(rateAMPua(:,1:85,:),[],2).*3));
            sig=(meancatdata)>thresh;
            sig=sig&validchannels;
            
        end
        
        
        for chn=chnrange
            for stimchn=1:size(sig,3)
                if sig(chn,:,stimchn)==1
                    alldata{ampit}=[alldata{ampit}; rateAMPua(chn,:,stimchn)];
                end
            end
        end
    end
end
end
