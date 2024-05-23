function weirdartifactsp=removeArtifactSpikes(sp)

    % Find minimum value and its index
    MinY = min(sp(:,2:end),[],2)/2;
    
    % Determine half Minimum value
    belowhalfMin = sp(:,2:end)<=MinY;
    
% Find first and last index where each row exceeds or equals half-min value
% from the minimum point at 13
first_indices = arrayfun(@(row) find([true ~belowhalfMin(row, 1:13)], 1, 'last'), 1:size(belowhalfMin, 1),'UniformOutput',true);
last_indices = arrayfun(@(row) find([~belowhalfMin(row, 13:end) true], 1, 'first'), 1:size(belowhalfMin, 1))+11;

%remove weird artifacts
halfWidth1 = last_indices - first_indices;

%also remove any large wide positive spikes
% Find max value and its index
    [MaxY, Imax] = max(sp(:,2:end),[],2);
    MaxY(MaxY<0)=0;
Halfmax=MaxY./2;
   % Determine half Minimum value
    abovehalfMax = sp(:,2:end)>=Halfmax;
    
% Find first and last index where each row exceeds or equals half-min value
% from the minimum point at 13
first_indices = arrayfun(@(row) find([true ~abovehalfMax(row, 1:Imax(row))], 1, 'last'), 1:size(abovehalfMax, 1),'UniformOutput',true);
last_indices = arrayfun(@(row) find([~abovehalfMax(row, Imax(row):end) true], 1, 'first'), 1:size(abovehalfMax, 1))+Imax'-1;

%remove weird artifacts
halfWidth2 = last_indices - first_indices;

%catch ends that are bad and indicative of artifact
[r,~]=find(sp(:,40:end)>50);

weirdartifactsp=[sp(r,:); sp(halfWidth1>10,:); sp(halfWidth2>10,:)];
weirdartifactsp=unique(weirdartifactsp(:,1));

end