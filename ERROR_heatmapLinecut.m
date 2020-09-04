function ERROR_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse)
%Plots true data and prediction, then plots error

TrueData_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,stderrspktrial,startpointseconds, secondstoanalyse)
maxid=originalEND/(n_REP_true*2);
endtrialelect=(find(cell2mat(trialinfo((trialjump*2+1):end,18))==-1,1)+trialjump*2-1)/2; %trials for one set of conditions 
if VarAmp==1
    number_trufig=maxid/(endtrialelect*length(NORECORDELECT));
else
    
end

AdditivePrediction_heatmapLinecut(AMPInterestSingleLinePlot,avgnospT,startpointseconds, secondstoanalyse)

end

