clc; 
clear;

% nsfile = '~/Documents/data/CJ221.noisegrid.044310.mat'; %CJ221
% nsfile = '~/Documents/data/cj223.noisegrid.115113.mat'; 
% nsfile = '/home/marmolab/data/2022/09/12/CJ223.noisegrid.213911.mat';
%nsfile = '/home/marmolab/data/2022/09/13/CJ223.noisegrid.145016.mat';
%nsfile = '/home/marmolab/data/2022/10/11/cj225.noisegrid.163819.mat';
nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\10\11\cj225.noisegrid.163819.mat';
%nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\10\11\cj225.noisegrid.113604.mat';
nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\11\22\cj230.noisegrid.154618.mat';
%nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\09\07\cj222.noisegrid.204533.mat';
%nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\10\11\cj225.noisegrid.103347.mat';
%nsfile='E:\DATA\MarmoblitzVisualDataCJ\2022\11\22\cj230.noisegrid.093513.mat';
%nsfile = 'E:\DATA\MarmoblitzVisualDataCJ\2022\09\07\cj222.noisegrid.115113.mat';
chList = 1:128;%[32,1,25,31,27,18,41,34,43,52,56,54,53,59,60,96,72,65,89,95,71,66,90,94,70,67,91,93,69,68,92];%1:128;%[128,13,76,41,88,23,5,120,89,103,82,99,16,70,56,38,105,49,79,15,40,1,37,80,45,53,18,122,58,115,61,121,74,64,104,25,106,117,62,11,31,78,90,98,20,68,125,46,17,52,59,2,6,93,8,29,26,81,102,77,54,12,107,109,47,84,55,19,111,7,108,91,97,34,27,42,60,69,57,75,83,100,63,10,113,87,51,92,9,33,36,123,44,118,127,116,119,85,35,86,32,95,67,24,101,126,4,73,94,65,3,124,50,110,72,43,71,22,28,112,21,114,96,48,66,30,14,39]; ; %[33   36   39    9   24   40   61  101   48   43   17   57  104   97   11  108   65   52   15   96  103  106   58   23   16  105   19  125   55   72   70  112  107   71   89   94  110   37   66   38   10  116   45   67   90   42   44   18]; %1:128;
fileS=dir;
alreadysaveddobject=find(contains({fileS.name},'dobject'));
if ~isempty(alreadysaveddobject)
    load(fileS(alreadysaveddobject).name)
else
    d = mapping.analysis.noisegrid(nsfile,'loadArgs', ...
        {'spikes',true,'source','ghetto','reload', true 'channels', chList});
end
    %%
    subtract_baseline=0; % if true uses grey squares for BL resp
totCh = length(chList); 
chArea = zeros(1,128);
chArea(chList<65) = 2;
chArea(chList>64) = 1;
% chArea = 1+ones(1, totCh);
% chArea(1) = 1; 
aSysOrder = randsample(1:totCh, totCh);
ROI = [0 15; -15 0];
aList = {'V1', 'V2'};

sWin = 60;%window to look in after randel change in ms
sBin = 50; %amount of smoothing for sdf - basically a movmean bin
offset=20;%ms offset from start of trial to start looking for spikes
interp = 1;
tickSp = 5;
figure(3); clf;
h = figure(4); clf;
RF = cell(1, totCh); 
firstTime = true; 
trialcutoff=0;
for ich = aSysOrder
    fprintf('Loading Channel %i ... \n', chList(ich)); 


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
    end
    
    [~, nSpikeTrials, ~] = size(d.spikes.spk); 
    
    nSpikesTot = sum(cellfun(@length, d.spikes.spk(:,:,ich))); %total number of spikes on that channel for all trials&frames
    fprintf('%i spikes detected ... \n', nSpikesTot); 
   
    allSpikes{ich} = zeros(1, nSpikesTot);
    spikeInd = ones(1, totCh); 
    
    fprintf('Getting spike times...')
    for iTrial = 1:nStimTrials
        nFrames = sum(nFramesList(:, iTrial));%total number of randel/square onsets for all randels in one trial
        
        if iTrial <= nSpikeTrials
            if iTrial>trialcutoff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%manually get rid of start trials that might be bad - NEED TO GET RID OF FIRST 25 TRIALS AND MAKE SURE THINGS ARE ALIGNED PROPERLY
              sTimes = round(d.spikes.spk{1,iTrial,ich}*1000); %spikeTimes in ms for 1 channel in 1 trial
              sCount = length(sTimes); 
              allSpikes{ich}(spikeInd(ich):spikeInd(ich)+sCount-1) = sTimes+tOffset; %spike time at start of trial plus 1ms
            else
                sTimes=[];
                sCount = length(d.spikes.spk{1,iTrial,ich}); 
            end
            spikeInd(ich) = spikeInd(ich) + sCount;%total number of spikes for each channel + number for this trial
        else
           warning('No spikes for this trial!');  
        end
  
%         if ich == 1
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
    allSpikes{ich} = allSpikes{ich}(allSpikes{ich}>0);%checks all the spikes were registered
    sTrain_sitmulus(allSpikes{ich}) = 1; %makes index of 1 for any randel onsets that had a spike 
     nSpikesTot= length(allSpikes{ich});
     if subtract_baseline==1
     %%%%%%%%%%%%Baseline calc
     tOffset = 1; 
    frameOffset = 0;
    if firstTime
        fprintf('Getting trial onsets...')
        [x_BL,y_BL,onsets] = d.getSingleSquareOnsets('greyonly',true); % onset times of the flashed squares
        
        %onsets [x*y,trials] - x*y is number randels or number squares
        %x= xvalues for whole grid. should be column 1 first i.e. all same
        %xvalue, then column 2 etc. 
        %y = yvalues for whole grid. should be column 1 first, i.e.
        %decreasing as you move down, then column 2 etc. 
        

        inROIx = x_BL >= ROI(1,1) & x_BL <= ROI(1,2);%determines if xval/yval within the ROI you declared above
        inROIy = y_BL >= ROI(2,1) & y_BL <= ROI(2,2);
        
        useRand_BL = inROIx & inROIy;%gets index of randels within the ROI
        xroi = x_BL(useRand_BL); yroi = y_BL(useRand_BL); 
        %onroi = onsets{inROIx & inROIy};

        onMS_baseline = cellfun(@(x) round(x.*1000), onsets, 'UniformOutput', false);%randel onset time in ms

        nFramesList_BL = cellfun(@length, onsets); %determines how many onsets in each trial for that randel
        nFramesList_BL = nFramesList_BL(useRand_BL, :); %only use the randels within the ROI
        nRand = sum(useRand_BL);
        [~, nStimTrials] = size(onsets);  
        
        onsetInds_baseline = zeros(3, sum(nFramesList_BL(:)));%total number of square onsets for all randels over all trials
        
    end

   
    
    fprintf('Getting spike times...')
    for iTrial = 1:nStimTrials
        nFrames = sum(nFramesList_BL(:, iTrial));%total number of randel/square onsets for all randels in one trial
        

  
%         if ich == 1
            timeList  = zeros(1, nFrames);
            idList    = zeros(1, nFrames);
            frameInd  = 1; 
            for iRand = find(useRand_BL)% iterate for randels in ROI
                t = onMS_baseline{iRand, iTrial}; %get one randel's onsets for current trial
                nf = length(t); %total number of onsets
                timeList(frameInd:frameInd+nf-1) = t;%saving those onsets in the array with total number of randel/square onsets for all randels in one trial
                idList(frameInd:frameInd+nf-1) = iRand*ones(1, nf);%which randel onsets have been saved in timelist already
                frameInd = frameInd + nf;%total number of onsets
            end   
    
            [orderTime, ind] = sort(timeList);%Timlist is currnently in order of randel location e.g. onsets for [1,1] first in the list then [2,1] etc. Need them in order of time. e.g. square 1 and square 20 might have flashed together first 
            orderId = idList(ind); %sort the randel list into time order too
            onsetInds_baseline(1, 1+frameOffset:nFrames+frameOffset) = orderTime + tOffset; %save the onset times of the randels
            onsetInds_baseline(2, 1+frameOffset:nFrames+frameOffset) = x_BL(orderId); % save the respective x coordinates for the randel list in order of onset time
            onsetInds_baseline(3, 1+ frameOffset:nFrames+frameOffset) = y_BL(orderId); % save the respective y coordinates for the randel list in order of onset time
            
            xs = unique(x_BL(useRand_BL)); %get only the unique x and y coordinates 
            ys = unique(y_BL(useRand_BL));
%         end
   
    
    tOffset = tOffset + max(timeList) + 500;%offset for the next trial 
    frameOffset = frameOffset + nFrames;%total number of onsets
    end
     else
         onsetInds_baseline=[];
    end
       
    
    

    %%%%%%%%%%%%RF plotting
%     RF{ich} = getRFMap(sTrain_sitmulus, onsetInds_stimulus,onsetInds_baseline, sWin, sBin, offset);
%    
    if nSpikesTot > 300%only include in the overlapping RF is there were more than 300 spikes
        RF{ich} = getRFMap(sTrain_sitmulus, onsetInds_stimulus,onsetInds_baseline, sWin, sBin, offset);
    else
        RF{ich} = zeros(length(ys),length(xs));
    end

        figure(7); 
        E_MAP=Depth(1);
        subplot(8,16, find(E_MAP==ich)); 
     
        imagesc(xs,ys,(RF{ich})); colormap('gray'); 
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        drawnow;
        set(gca,'YDir','normal')
    
        drawOverlappingRFs(h, RF, chArea, xs, ys, tickSp, aList);
    
    firstTime = false;
end


%% plot all RFs
E_MAP=Depth(1);
   arrayshape=reshape(E_MAP,16,8);
   arrayshape=arrayshape(:,[1 3 4 2 5 7 8 6]);
   arrayshape=flipud(arrayshape);
     for ich=1:128
        figure(8); 
        [r,c]=find(arrayshape==ich);
        subplot(16,8, (r - 1) * 8 + c); 
     
        imagesc(xs,ys,(RF{ich})); colormap('gray'); 
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        drawnow;
        set(gca,'YDir','normal')
     end
 
%%