function RF=getRFsinglechn(d,chn,timecourse,sWin)
chList=1:128;
totCh = length(chList); 
chArea = zeros(1,128);
chArea(chList<65) = 2;
chArea(chList>64) = 1;
% chArea = 1+ones(1, totCh);
% chArea(1) = 1; 
ROI = [0 15; -15 0];
aList = {'V1', 'V2'};

sBin = 50; %amount of smoothing for sdf - basically a movmean bin
interp = 1;
tickSp = 5;
RF = cell(1, length(timecourse)); 
firstTime = true; 
trialcutoff=0;
for timeit=1:length(timecourse)
    fprintf('Loading Channel %i ... \n', chn); 
    tOffset = 1; 
    frameOffset = 0;
    if firstTime
        fprintf('Getting trial onsets...')
        [x,y,onsets] = d.getSingleSquareOnsets(); % onset times of the flashed squares
        
        %onsets [x*y,trials] - x*y is number randels or number squares
        %x= xvalues for whole grid. should be column 1 first i.e. all same
        %xvalue, then column 2 etc. 
        %y = yvalues for whole grid. should be column 1 first, i.e.
        %decreasing as you move down, then column 2 etc. 
        

        inROIx = x >= ROI(1,1) & x <= ROI(1,2);%determines if xval/yval within the ROI you declared above
        inROIy = y >= ROI(2,1) & y <= ROI(2,2);
        
        useRand = inROIx & inROIy;%gets index of randels within the ROI
        xroi = x(useRand); yroi = y(useRand); 
        %onroi = onsets{inROIx & inROIy};

        onMS = cellfun(@(x) round(x.*1000), onsets, 'UniformOutput', false);%randel onset time in ms

        nFramesList = cellfun(@length, onsets); %determines how many onsets in each trial for that randel
        nFramesList = nFramesList(useRand, :); %only use the randels within the ROI
        nRand = sum(useRand);
        [~, nStimTrials] = size(onsets);  
        
        onsetInds_stimulus = zeros(3, sum(nFramesList(:)));%total number of square onsets for all randels over all trials
        
        allSpikes = cell(1, d.spikes.numChannels); 
    
    
    [~, nSpikeTrials, ~] = size(d.spikes.spk); 
    
    nSpikesTot = sum(cellfun(@length, d.spikes.spk(:,:,chn))); %total number of spikes on that channel for all trials&frames
    fprintf('%i spikes detected ... \n', nSpikesTot); 
   
    allSpikes{chn} = zeros(1, nSpikesTot);
    spikeInd = ones(1, totCh); 
    
    fprintf('Getting spike times...')
    for iTrial = 1:nStimTrials
        nFrames = sum(nFramesList(:, iTrial));%total number of randel/square onsets for all randels in one trial
        
        if iTrial <= nSpikeTrials
            if iTrial>trialcutoff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%manually get rid of start trials that might be bad - NEED TO GET RID OF FIRST 25 TRIALS AND MAKE SURE THINGS ARE ALIGNED PROPERLY
              sTimes = round(d.spikes.spk{1,iTrial,chn}*1000); %spikeTimes in ms for 1 channel in 1 trial
              sCount = length(sTimes); 
              allSpikes{chn}(spikeInd(chn):spikeInd(chn)+sCount-1) = sTimes+tOffset; %spike time at start of trial plus 1ms
            else
                sTimes=[];
                sCount = length(d.spikes.spk{1,iTrial,chn}); 
            end
            spikeInd(chn) = spikeInd(chn) + sCount;%total number of spikes for each channel + number for this trial
        else
           warning('No spikes for this trial!');  
        end
  
%         if chn == 1
            timeList  = zeros(1, nFrames);
            idList    = zeros(1, nFrames);
            frameInd  = 1; 
            for iRand = find(useRand)% iterate for randels in ROI
                t = onMS{iRand, iTrial}; %get one randel's onsets for current trial
                nf = length(t); %total number of onsets
                timeList(frameInd:frameInd+nf-1) = t;%saving those onsets in the array with total number of randel/square onsets for all randels in one trial
                idList(frameInd:frameInd+nf-1) = iRand*ones(1, nf);%which randel onsets have been saved in timelist already
                frameInd = frameInd + nf;%total number of onsets
            end   
    
            [orderTime, ind] = sort(timeList);%Timlist is currnently in order of randel location e.g. onsets for [1,1] first in the list then [2,1] etc. Need them in order of time. e.g. square 1 and square 20 might have flashed together first 
            orderId = idList(ind); %sort the randel list into time order too
            onsetInds_stimulus(1, 1+frameOffset:nFrames+frameOffset) = orderTime + tOffset; %save the onset times of the randels
            onsetInds_stimulus(2, 1+frameOffset:nFrames+frameOffset) = x(orderId); % save the respective x coordinates for the randel list in order of onset time
            onsetInds_stimulus(3, 1+ frameOffset:nFrames+frameOffset) = y(orderId); % save the respective y coordinates for the randel list in order of onset time
            
            xs = unique(x(useRand)); %get only the unique x and y coordinates 
            ys = unique(y(useRand));
%         end
   
    
    tOffset = tOffset + max(timeList) + 500;%offset for the next trial 
    frameOffset = frameOffset + nFrames;%total number of onsets
   
    end
    
    sTrain_sitmulus = zeros(1, tOffset+1000); 
    allSpikes{chn} = allSpikes{chn}(allSpikes{chn}>0);%checks all the spikes were registered
    sTrain_sitmulus(allSpikes{chn}) = 1; %makes index of 1 for any randel onsets that had a spike 
     nSpikesTot= length(allSpikes{chn});
     onsetInds_baseline=[];
       
    
    end

    %%%%%%%%%%%%RF plotting
%     RF{chn} = getRFMap(sTrain_sitmulus, onsetInds_stimulus,onsetInds_baseline, sWin, sBin, offset);
%    
    if nSpikesTot > 300%only include in the overlapping RF is there were more than 300 spikes
        RF{timeit} = getRFMap(sTrain_sitmulus, onsetInds_stimulus,onsetInds_baseline, sWin, sBin, timecourse(timeit));
    else
        RF{timeit} = zeros(length(ys),length(xs));
    end
%     
        figure(10); 
        hold on
        E_MAP=Depth(1);
        subplot(2,ceil(length(timecourse)/2), timeit); 
     
        imagesc(xs,ys,(RF{timeit})); colormap('gray'); 
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        drawnow;
    if firstTime
        allRFs=zeros(length(ys),length(xs));
    end
    allRFs=allRFs+RF{timeit};
        %drawOverlappingRFs(h, RF, chArea, xs, ys, tickSp, aList);
    
    firstTime = false;
end
figure
imagesc(xs,ys,(allRFs)); colormap('gray');
end
