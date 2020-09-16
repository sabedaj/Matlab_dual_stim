load('overallsigstim.mat')
avgnum=numsigchn./stimelectrodes;

elect=find(~isnan(avgnum));
maxelect=find(avgnum==nanmax(avgnum));
endfrompeak=avgnum(maxelect:end);
beginfrompeak=[NaN; flipud(avgnum(1:maxelect-1))];
if length(endfrompeak)>length(beginfrompeak)
    diff=length(endfrompeak)-length(beginfrompeak);
    beginfrompeak=[beginfrompeak; NaN(diff,1)];
elseif length(endfrompeak)<length(beginfrompeak)
    diff=length(beginfrompeak)-length(endfrompeak);
    endfrompeak=[endfrompeak; NaN(diff,1)];
end
distfrompeak=nansum([beginfrompeak endfrompeak],2);
equal=find((~isnan(beginfrompeak)==~isnan(endfrompeak)));
equal=intersect(find(~isnan(endfrompeak)==1),equal);
distfrompeak(equal,:)=distfrompeak(equal,:)./2;
avgnum(isnan(avgnum))=[];
figure
plot(elect,avgnum)

figure
plot([0 50.*(1:(length(distfrompeak)-1))], distfrompeak)
title(['Distance from electrode ' num2str(maxelect) ' with peak activation'])

