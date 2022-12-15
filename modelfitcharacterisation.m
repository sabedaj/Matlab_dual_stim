function [rsq,AIC,MSE]=modelfitcharacterisation(y,yfit,yinput,numparam)
%y realdata
%yinput inputdata1
%yfit modeldata
% numparam number of input parameters to the model
NumObservations=sum(~isnan(yinput),'all');
DFE=NumObservations-numparam;%degrees freedom
MSE=(sum((y-yfit).^2,'all','omitnan'))/sum(~isnan(y),'all');%mean square error
Var = DFE/NumObservations*MSE;
w_r = ones(NumObservations,1);
L = -(DFE + NumObservations*log(2*pi) + sum(log(Var./w_r)))/2;
%L = -( sum((y-yfit).^2,'all','omitnan')/Var + NumObservations*log(2*pi) + NumObservations*log(Var))/2;%-ve loglikelihood
AIC=2*(numparam-L);
rsq=1-(sum((y-yfit).^2,'all','omitnan')./sum((y-nanmean(y,'all')).^2,'all','omitnan')); %need to check this https://au.mathworks.com/help/stats/coefficient-of-determination-r-squared.html
end