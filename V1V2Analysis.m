function alldata=V1V2Analysis(ratestruct,savefilename,chnrange,excitesupress, chnexclude)
%refactor epochratespk
%loop through data recorded at multiple current levels
%plot V1 and V2 responses

%ratestruct is a structure with fields that each represent a folder
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
AMPchosen=[2 5 6 8 10];

%subselect baseline
subselectbaseline=1:89;
permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
subselectbaseline=permuted_numbers(1:50);
ratespiking=cell(numfolderstotal,1);
alldata=cell(length(AMPchosen),1);
alldata_stim=cell(length(AMPchosen),1);


for ampit=1:length(AMPchosen)
    %loop through folders
    for folderit=1:numfolderstotal
        
        %loop through trial repeats and average
        trialend=length(fieldnames(ratestruct{numerfolders}));
        for trial=1:trialend
            trialcheck=['T' num2str(trial)];
            ratespiking{numerfolders}(:,:,trial)=nanmean(ratestruct{numerfolders}.(trialcheck),3);
        end
        %find the stimualtion with the nearest current to AMPchosen
        [~,i]=min(abs(savefilename{numerfolders}{4}(:,2)-ampinterest)); ampinterest=savefilename{numerfolders}{4}(i,2);%nearest neighbour amp
        
        rateAMPua=ratespiking{numerfolders}(:,:,savefilename{numerfolders}{4}((savefilename{numerfolders}{4}(:,2)==ampinterest),1));
        rateAMPua(chnexclude,:,:)=nan;
        %determine if we are evaluating suppression or excitation
        %excitation = 1 suppression =0
        %determine which channels are significantly responding
        if excitesupress==1
            %excitation
            %find channels that are significantly responding
            thresh=squeeze(mean(rateAMPua(:,subselectbaseline,:),2)+(std(rateAMPua(:,subselectbaseline,:),[],2).*3));
            sigchn=squeeze(mean(rateAMPua(:,92:181,:),2))>thresh;
        else
            %suppression
            %find channels that are significantly responding
            thresh=squeeze(mean(rateAMPua(:,subselectbaseline,:),2));
            sigchn=squeeze(mean(rateAMPua(:,90+25:181,:),2))<thresh;
        end
        %loop thorugh channels and stim channels
        for chn=chnrange
            for stimchn=1:size(sigchn,3)
                %pool significant channels from rateAMPua
                if sigchn(chn,stimchn)
                    %pool data from chosen chnrange
                    alldata{ampit}(iteratealldat,1:181)=rateAMPua(chn,:,stimchn);
               iteratealldat=iteratealldat+1;
               %pool stim data
               alldata_stim{ampit}(iteratestimdat:iteratestimdat+63,1:181)=rateAMPua(65:128,:,stimchn);
               iteratestimdat=iteratestimdat+64;

                end
            end
        end
        
        
        
    end
end





end