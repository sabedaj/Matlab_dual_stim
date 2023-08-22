function multielectplot(ratestruct,savefilename,chnrange,excitesupress, chnexclude)
ampinterest=6;
 s = RandStream('mt19937ar','Seed',296);
% pre-processing
numfolderstotal=size(ratestruct,1)-sum(cellfun(@isempty, ratestruct));
ratespiking=cell(numfolderstotal,1);

for numerfolders=1:numfolderstotal
    trialend=size(savefilename{numerfolders}{2}.AMP,1);
    for trial=1:trialend
        trialcheck=['T', num2str(trial)];
        ratespiking{numerfolders}(:,:,trial)=nanmean(ratestruct{numerfolders}.(trialcheck),3);
        %make sure only responding elect included
        subselectbaseline=1:89;
        permuted_numbers = subselectbaseline(randperm(s,length(subselectbaseline)));
        subselectbaseline=permuted_numbers(1:50);
        thresh=squeeze(mean(ratespiking{numerfolders}(:,subselectbaseline,trial),2)+(std(ratespiking{numerfolders}(:,subselectbaseline,trial),[],2).*3));
        notsigsingle=squeeze(mean(ratespiking{numerfolders}(:,92:100,trial),2))<=thresh;
        ratespiking{numerfolders}(notsigsingle,:,trial)=nan;
    end
end

%%
% 1. # elect vs FR.
numelecttotal=length(savefilename{1}{2}.CHN);
separateeachstimelect=cell(numelecttotal,numelecttotal);% 8 plots with 8 different #elect
respelect=cell(numelecttotal,1);
V2resp=cell(numelecttotal,1);
V1resp=cell(numelecttotal,1);
for numstimelect=1:numelecttotal
    for numerfolders=1:numfolderstotal
        trialend=size(savefilename{numerfolders}{2}.AMP,1);
        for trial=1:trialend
            whichelectstim=savefilename{numerfolders}{2}.AMP(trial,:)==ampinterest;
            if sum(whichelectstim)~=numstimelect%check if the trial has the correct num stim elect and correct amplitude stim
                continue
            end
            index=find(whichelectstim);
          for i=1:numstimelect
            separateeachstimelect{index(i),numstimelect}=cat(1,separateeachstimelect{index(i),numstimelect},ratespiking{numerfolders}(1:64,:,trial));
          end
            diff(squeeze(mean(ratespiking{numerfolders}(:,92:181,:),2)));
            V2resp{numstimelect}=cat(1,V2resp{numstimelect},ratespiking{numerfolders}(1:64,:,trial));
            V1resp{numstimelect}=cat(1,V1resp{numstimelect},ratespiking{numerfolders}(65:128,:,trial));
            respelect{numstimelect}=cat(3,respelect{numstimelect},ratespiking{numerfolders}(:,:,trial));
        end
    end
end
figure
ax=axes;
hold on
color1 = linspace(0,1,numelecttotal);
newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
legendCell=cell(numelecttotal,1);
for numelect=1:numelecttotal
   stdshade(V1resp{numelect}.*1000,0.2,newcolors(numelect,:),[-90:90],1,ax); 
   legendCell{numelect} = num2str(numelect);
end
xlim([-50 85])
leg=legend(legendCell);
title(leg,'# stim elect')
xlabel('Time (ms)')
ylabel('Firing rate (sp/s)')
title('V1')
beautifyPlot;


figure
ax=axes;
hold on
color1 = linspace(0,1,numelecttotal);
newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
legendCell=cell(numelecttotal,1);
for numelect=1:2:numelecttotal
   stdshade(V2resp{numelect}.*1000,0.2,newcolors(numelect,:),[-90:90],1,ax); 
   legendCell{numelect} = num2str(numelect);
end
xlim([-50 85])
leg=legend(legendCell);
title(leg,'# stim elect')
xlabel('Time (ms)')
ylabel('Firing rate (sp/s)')
title('V2')
beautifyPlot;

%%
for whichelectfocus=1:numelecttotal
figure
ax=axes;
hold on
color1 = linspace(0,1,numelecttotal);
newcolors = [flipud(color1') zeros(length(color1),1) (color1')];
legendCell=cell(numelecttotal,1);
for numelect=1:numelecttotal
   stdshade(separateeachstimelect{whichelectfocus,numelect}.*1000,0.2,newcolors(numelect,:),[-90:90],1,ax); 
   legendCell{numelect} = num2str(numelect);
end
xlim([-50 85])
leg=legend(legendCell);
title(leg,'# stim elect')
xlabel('Time (ms)') 
ylabel('Firing rate (sp/s)')
title(['V2 elect ' num2str(whichelectfocus)])
beautifyPlot;
end
%%
% 2. # elect vs spread

% 3. # elect vs distribution early late, large small



end