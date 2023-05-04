function [RF] = getRFMap(sTrain, onsetInds, sWin, sBin)

xs = unique(onsetInds(2,:)); nx = length(xs);
ys = unique(onsetInds(3,:)); ny = length(ys);
 
sdf = conv(sTrain, ones(1, sBin));
thisSDF = cell(length(xs), length(ys));
nTimes = length(sdf); 

% collect the SDFs for each channel/position
for x = 1:length(xs)
    for y = 1:length(ys)
        theseTimes = onsetInds(1, onsetInds(2,:) == xs(x) &...
                                  onsetInds(3,:) == ys(y)); 
        thisSDF{x,y} = zeros(length(theseTimes), sWin);
        for t = 1:length(theseTimes)
            if theseTimes(t) + sWin-1 < nTimes
                thisSDF{x,y}(t,:) = sdf(theseTimes(t):theseTimes(t) + sWin-1); 
            end
        end
    end
end
   
    % average across trials
    meanSDF = cellfun(@(x) squeeze(mean(x)), thisSDF, 'UniformOutput', false);
    
    % take the peak of the SDF
    RF = cellfun(@max, meanSDF);