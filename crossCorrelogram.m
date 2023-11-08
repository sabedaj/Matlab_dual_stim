function [ccg, lags] = crossCorrelogram(spikeTrain1, spikeTrain2, windowMs)
   % Compute the cross-correlogram of two spike trains in milliseconds
    
    % Input:
    % - spikeTrain1: A vector of spike times for neuron 1 (sorted) in milliseconds
    % - spikeTrain2: A vector of spike times for neuron 2 (sorted) in milliseconds
    % - windowMs: The size of the window to consider (in milliseconds)


    % Output:
    % - xcorr: A vector of the cross-correlation values
    % - lags: A vector of the lag values (in milliseconds) for each of the cross-correlation values
    ccg=[];
    lags=[];
    %iterate over first spike train
    for spikenum= 1 : length(spikeTrain1)
       % Calculate the time differences between spikes in the two spike trains
        timeDiffs = spikeTrain2 - spikeTrain1(spikenum);
        
        % Filter time differences within the specified range
        validDiffs = timeDiffs(timeDiffs >= -windowMs & timeDiffs <= windowMs);

        %hist in 1 ms bins
        [counts, lags] = histcounts(validDiffs, -windowMs:1:windowMs);
        %store counts and lags in matrix
        if spikenum == 1
            ccg = counts;
        else
            ccg = ccg + counts;
        end
    end
   
    %ccg = ccg / (length(spikeTrain1)+length(spikeTrain2));%normalize by number of spikes in train1

end
