function [RF] = getRFMap(sTrain, onsetInds,onsetInds_baseline, sWin, sBin, timeoffset)
%spike counter
xs = unique(onsetInds(2,:)); nx = length(xs);
ys = unique(onsetInds(3,:)); ny = length(ys);
 
sdf = conv(sTrain, ones(1, sBin));%smooth out spike train
thisSDF = cell(length(ys), length(xs));
thisSDF_BL = cell(length(ys), length(xs));
nTimes = length(sdf); 

% collect the SDFs for each channel/position
for x = 1:length(xs)
    for y = 1:length(ys)
        theseTimes = onsetInds(1, onsetInds(2,:) == xs(x) &...
            onsetInds(3,:) == ys(y)); %get onset times for a specific randel
        
        
        thisSDF{y,x} = zeros(length(theseTimes), sWin);
        for t = 1:length(theseTimes)
            if theseTimes(t) + sWin-1 < nTimes
                thisSDF{y,x}(t,:) = sdf(theseTimes(t) + timeoffset:theseTimes(t)+ timeoffset + sWin-1); % for each onset with time offset, get the spikes in the chosen window
                
            end
        end
        if ~isempty(onsetInds_baseline)
            theseTimes_BL = onsetInds_baseline(1, onsetInds_baseline(2,:) == xs(x) &...
                onsetInds_baseline(3,:) == ys(y)); %get onset times for a specific randel
            thisSDF_BL{y,x} = zeros(length(theseTimes_BL), sWin);
            for t = 1:length(theseTimes_BL)
                if theseTimes_BL(t) + sWin-1 < nTimes
                    thisSDF_BL{y,x}(t,:) = sdf(theseTimes_BL(t) + timeoffset:theseTimes_BL(t)+ timeoffset + sWin-1); % for each onset with time offset, get the spikes in the chosen window
                    
                end
            end
        end
        
        
    end
end
    
    % average across trials
    if ~isempty(onsetInds_baseline)
    meanSDF_sig = cellfun(@(x) squeeze(mean(x)), thisSDF, 'UniformOutput', false);%for every randel onset, we have the chosen window (ms) of data. Average across randel onsets to get one signal for the chosen window e.g. if window is 100ms, 100 elements representing 1ms each will be present
    meanSDF_BL= cellfun(@(x) squeeze(mean(x)), thisSDF_BL, 'UniformOutput', false);
    meanSDF = cellfun(@(x,y) x-y,meanSDF_sig,meanSDF_BL, 'UniformOutput', false);
    else
        meanSDF  = cellfun(@(x) squeeze(mean(x)), thisSDF, 'UniformOutput', false);%for every randel onset, we have the chosen window (ms) of data. Average across randel onsets to get one signal for the chosen window e.g. if window is 100ms, 100 elements representing 1ms each will be present
    end
    % take the peak of the SDF
    RF = cellfun(@max, meanSDF);% get the max firing rate in the window
   %  RF = cellfun(@mean, meanSDF);% get the max firing rate in the window