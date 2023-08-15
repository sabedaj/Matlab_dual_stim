clc; 
clear;

% nsfile = '~/Documents/data/CJ221.noisegrid.044310.mat'; %CJ221
% nsfile = '~/Documents/data/cj223.noisegrid.115113.mat'; 
% nsfile = '/home/marmolab/data/2022/09/12/CJ223.noisegrid.213911.mat';
%nsfile = '/home/marmolab/data/2022/09/13/CJ223.noisegrid.145016.mat';
%nsfile = '/home/marmolab/data/2022/10/11/cj225.noisegrid.163819.mat';
%nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\10\11\cj225.noisegrid.163819.mat';
nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\09\07\cj222.noisegrid.115113.mat';
chList = [128,13,76,41,88,23,5,120,89,103,82,99,16,70,56,38,105,49,79,15,40,1,37,80,45,53,18,122,58,115,61,121,74,64,104,25,106,117,62,11,31,78,90,98,20,68,125,46,17,52,59,2,6,93,8,29,26,81,102,77,54,12,107,109,47,84,55,19,111,7,108,91,97,34,27,42,60,69,57,75,83,100,63,10,113,87,51,92,9,33,36,123,44,118,127,116,119,85,35,86,32,95,67,24,101,126,4,73,94,65,3,124,50,110,72,43,71,22,28,112,21,114,96,48,66,30,14,39]; ; %[33   36   39    9   24   40   61  101   48   43   17   57  104   97   11  108   65   52   15   96  103  106   58   23   16  105   19  125   55   72   70  112  107   71   89   94  110   37   66   38   10  116   45   67   90   42   44   18]; %1:128;
totCh = length(chList); 
chArea = zeros(1,128);
chArea(chList<65) = 2;
chArea(chList>64) = 1;
% chArea = 1+ones(1, totCh);
% chArea(1) = 1; 
aSysOrder = randsample(1:totCh, totCh);
ROI = [-10 0; 0 10];
aList = {'V1', 'V2'};

sWin = 100;
sBin = 50; 
interp = 1;
tickSp = 5;
figure(1); clf;
h = figure(2); clf;
RF = cell(1, totCh); 
firstTime = true; 
for ich = aSysOrder
    fprintf('Loading Channel %i ... \n', chList(ich)); 
    d = mapping.analysis.noisegrid(nsfile,'loadArgs', ...
        {'spikes',true,'source','ghetto','reload', false, 'channels', chList(ich)});%source ghetto loads raw data every time
       % d=mapping.analysis.noisegrid('cj222.noisegrid.115113.mat',{'loadEye',false,'spikes',true,'source','ghetto','reload',true,'cfg','acute.H64Flexi64FlexiIntan','channels', chList(ich),'useCAR',false});
    
    
    tOffset = 1; 
    frameOffset = 0;
    if firstTime
        fprintf('Getting trial onsets...')
        [x,y,onsets] = d.getSingleSquareOnsets(); % onset times of the flashed squares
        

        inROIx = x >= ROI(1,1) & x <= ROI(1,2);
        inROIy = y >= ROI(2,1) & y <= ROI(2,2);
        
        useRand = inROIx & inROIy;
        xroi = x(useRand); yroi = y(useRand); 
        %onroi = onsets{inROIx & inROIy};

        onMS = cellfun(@(x) round(x.*1000), onsets, 'UniformOutput', false);

        nFramesList = cellfun(@length, onsets);
        nFramesList = nFramesList(useRand, :);
        nRand = sum(useRand);
        [~, nStimTrials] = size(onsets);  
        
        onsetInds = zeros(3, sum(nFramesList(:)));
        
        allSpikes = cell(1, d.spikes.numChannels); 
    end
    
    [~, nSpikeTrials, ~] = size(d.spikes.spk); 
    
    nSpikesTot = sum(cellfun(@length, d.spikes.spk)); 
    fprintf('%i spikes detected ... \n', nSpikesTot); 
   
    allSpikes{ich} = zeros(1, nSpikesTot);
    spikeInd = ones(1, totCh); 
    
    fprintf('Getting spike times...')
    for iTrial = 1:nStimTrials
        nFrames = sum(nFramesList(:, iTrial));
        
        if iTrial <= nSpikeTrials
            sTimes = round(d.spikes.spk{iTrial}*1000); %spikeTimes in ms 
            sCount = length(sTimes); 
            allSpikes{ich}(spikeInd(ich):spikeInd(ich)+sCount-1) = sTimes+tOffset; 
            spikeInd(ich) = spikeInd(ich) + sCount;
        else
           warning('No spikes for this trial!');  
        end
  
%         if ich == 1
            timeList  = zeros(1, nFrames); 
            idList    = zeros(1, nFrames);
            frameInd  = 1; 
            for iRand = find(useRand)
                t = onMS{iRand, iTrial}; 
                nf = length(t); 
                timeList(frameInd:frameInd+nf-1) = t;
                idList(frameInd:frameInd+nf-1) = iRand*ones(1, nf);
                frameInd = frameInd + nf;
            end   
    
            [orderTime, ind] = sort(timeList);
            orderId = idList(ind); 
            onsetInds(1, 1+frameOffset:nFrames+frameOffset) = orderTime + tOffset;
            onsetInds(2, 1+frameOffset:nFrames+frameOffset) = x(orderId); 
            onsetInds(3, 1+ frameOffset:nFrames+frameOffset) = y(orderId);
            
            xs = unique(x(useRand)); 
            ys = unique(y(useRand));
%         end
   
    
    tOffset = tOffset + max(timeList) + 500; 
    frameOffset = frameOffset + nFrames;
   
    end
    
    sTrain = zeros(1, tOffset+1000); 
    allSpikes{ich} = allSpikes{ich}(allSpikes{ich}>0);
    sTrain(allSpikes{ich}) = 1; 
    if nSpikesTot > 300
        RF{ich} = getRFMap(sTrain, onsetInds, sWin, sBin);
    else
        RF{ich} = zeros(length(xs), length(ys));
    end
    
        figure(1); 
        subplot(11,12, ich); 
     
        imagesc(flipud(RF{ich})); colormap('gray'); 
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        drawnow;
    
    
        drawOverlappingRFs(h, RF, chArea, xs, ys, tickSp, aList);
    
    firstTime = false;
end


%%