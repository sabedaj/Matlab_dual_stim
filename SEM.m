function [standarderrormean]=SEM(data,rowflag)
%used to calculate the standard error of input data. Regardless of array or
%vector it does all elements
if rowflag==0
standarderrormean=std(data,0,'all','omitnan')/sqrt(sum(~isnan(data),'all'));
else
    standarderrormean=std(data,0,2,'omitnan')/sqrt(size(data,2));
end
end