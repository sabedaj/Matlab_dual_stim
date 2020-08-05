function medianTrace = applyCARtoDat(filename, nChansTotal, outputDir, dName)
% Subtracts median of each channel, then subtracts median of each time
% point.
%
% filename should include the extension
% outputDir is optional, by default will write to the directory of the input file
%
% should make chunk size as big as possible so that the medians of the
% channels differ little from chunk to chunk.
chunkSize = 3e4*128;

fid = []; fidOut = [];

d = dir(filename);
nSampsTotal = d.bytes/nChansTotal/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);
try
  
  [pathstr, name, ext] = fileparts(filename);
  fid = fopen(filename, 'r');
  if nargin < 3
    outputFilename  = [pathstr filesep name '_CAR' ext];
  else
    outputFilename  = [outputDir filesep name '_CAR' ext];
  end
  fidOut = fopen(outputFilename, 'w');
  
  % theseInds = 0;
  chunkInd = 1;
  medianTrace = zeros(1, nSampsTotal);
  dispstat('','init');
  dispstat(sprintf('Processing CAR . . .'),'keepthis','n');    
  while 1
    dispstat(sprintf('Progress %03.2f%%',(100*(chunkInd/nChunksTotal))),'timestamp');
    if strcmp(dName,'analogin')
        dat = fread(fid, [nChansTotal chunkSize], 'uint16');
        dat = (dat - 32768) .* 0.0003125;
    elseif strcmp(dName,'amplifier')
        dat = fread(fid, [nChansTotal chunkSize], 'int16') .* 0.195;
    end
    
    if ~isempty(dat)
      
      %         theseInds = theseInds(end):theseInds(end)+chunkSize-1;
      
      dat = bsxfun(@minus, dat, median(dat,2)); % subtract median of each channel
      tm = median(dat,1);
      dat = bsxfun(@minus, dat, tm); % subtract median of each time point
      dat = dat ./ 0.195;
      fwrite(fidOut, dat, 'int16');
      medianTrace((chunkInd-1)*chunkSize+1:(chunkInd-1)*chunkSize+numel(tm)) = tm;
      
    else
      break
    end
    
    chunkInd = chunkInd+1;
  end
  
  fclose(fid);
  fclose(fidOut);
  
catch me
  
  if ~isempty(fid)
    fclose(fid);
  end
  
  if ~isempty(fidOut)
    fclose(fidOut);
  end
  
  
  rethrow(me)
  
end
