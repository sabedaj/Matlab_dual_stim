e=25;%maximum number electrodes included for white noise stim e.g. 8 for up to 8 electrodes
current=1; %number levels e.g. 1:2 for 2 current levels
repeats=1; % number of repeats
n=1:e;
row=0;

array=zeros(n(end)*2+1,n(end));
for k=1:e
    permutes=nchoosek(n,k);
    allcurrents = repelem(current, size(permutes,2));%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:size(permutes,1)
        currentcombo=uniqueperms(allcurrents);
        currentcombo=unique(currentcombo(:,1:k),'rows');
        for currentcycle=1:size(currentcombo,1)
            row=row+1;
            array(row,permutes(i,:))=currentcombo(currentcycle,:);
        end
    end
end
array(row+1,:)=0;
fprintf(['Exp will take ' num2str((size(array,1)*repeats)/3600) ' hours\n']) %time in hours