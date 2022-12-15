function smoothedenvelope=smoothData(data, FilterLength)
%smoothing parameters
b=ones(FilterLength,1)./FilterLength;
smoothedenvelope=nan(32,size(data,2));
for paircount=1:size(data,2)
    datnonan=data(:,paircount);
    endnonan=find(diff(isnan(datnonan))==1);
    startnonan=find(diff(isnan(datnonan))==-1);
    if isempty(startnonan)
        startnonan=0;
    end
    [~,p]=max(endnonan-startnonan);
    lastnonan=endnonan(p);
    firstnonan=startnonan(p)+1;
    datnonan=datnonan(firstnonan:lastnonan);
    %                     if length(datnonan)<7
    %                         continue
    %                     end
    if startnonan==0 && isempty(datnonan)
        continue
    end
    smoothedenvelopesingle=filtfilt(b,1,datnonan);
    smoothedenvelope(:,paircount)=[nan(firstnonan-1,1); smoothedenvelopesingle; nan(32-lastnonan,1)];
end
end