function [centroidY,centroidX]=centroidpos(data)
%calculates the centroid of data array



data(isnan(data))=0;
data(data<0)=0;
% Get the size of the array
[rows, cols] = size(data);

% Initialize variables to store the weighted sums
weightedSumX = 0;
weightedSumY = 0;
totalWeight = 0;

% Loop through all elements of the array
for row = 1:rows
    for col = 1:cols
        % Get the value at the current position
        value = data(row, col);
        
        % Calculate the weighted sum of x and y coordinates
        weightedSumX = weightedSumX + col * value;
        weightedSumY = weightedSumY + row * value;
        
        % Update the total weight
        totalWeight = totalWeight + value;
    end
end

% Calculate the centroid coordinates
centroidX = weightedSumX / totalWeight;
centroidY = weightedSumY / totalWeight;




end