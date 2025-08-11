function [directions, directionCounts, directionLabels, clusterID, eventCounts, eventCounts_clusterID] = clusterWaves_eventCount6_CSD(data, cellTypes)
% classifyCSDResponse Classifies CSD responses into chronic and acute clusters
%
% Input:
%   data - Nx3 matrix where:
%          Column 2 = preCSD
%          Column 3 = duringCSD
%          Column 4 = postCSD
%
% Outputs:
%   chronicCluster - Nx1 vector with values:
%                    1  = increase chronic
%                   -1  = decrease chronic
%                    0  = no change chronic
%
%   acuteCluster   - Nx1 vector with values:
%                    1 = increase acute
%                    0 = no change acute
%
%   combinedLabel  - Nx1 categorical array describing both classifications
%                    (e.g., "Chronic↑ + Acute↑")


pre = data(:, 2);
during = data(:, 3);
post = data(:, 4);

eventCounts = [pre, during, post];

% PRE-TO-DURING Logical conditions for classification of events 
increasingEventsChange1 = eventCounts(:, 2) >= 1;
noChangeEventsChange1 = ~(increasingEventsChange1);

% PRE-TO-POST Logical conditions for classification of events 
increasingEventsChange2 = eventCounts(:, 3) >= eventCounts(:, 1) + 2;
decreasingEventsChange2 = eventCounts(:, 3) <= eventCounts(:, 1) - 2;
noChangeEventsChange2 = ~(increasingEventsChange2 | decreasingEventsChange2);

labels1 = repmat("no change", size(eventCounts, 1), 1);
labels1(increasingEventsChange1) = "increasing";
labels1(noChangeEventsChange1) = "noChange";

labels2 = repmat("no change", size(eventCounts, 1), 1);
labels2(increasingEventsChange2) = "increasing";
labels2(decreasingEventsChange2) = "decreasing";

% Preallocate
directionLabels = strings(size(eventCounts,1), 1)';

% Create 9 clusters based on change1 and change2
for i = 1:length(eventCounts)
    acute = labels1(i);  % pre → during
    chronic = labels2(i);    % pre → post

    if acute == "increasing"
        if chronic == "increasing"
            directionLabels(i) = "↑↑";
        elseif chronic == "no change"
            directionLabels(i) = "↑x";
        elseif chronic == "decreasing"
            directionLabels(i) = "↑↓";
        end
    elseif acute == "noChange"
        if chronic == "increasing"
            directionLabels(i) = "x↑";
        elseif chronic == "no change"
            directionLabels(i) = "xx";
        elseif chronic == "decreasing"
            directionLabels(i) = "x↓";
        end
    end
end

directions = ["↑↑", "↑x", "↑↓", "x↑", "xx", "x↓"];
numDirections = numel(directions);

% Unique cell types
[uniqueTypes, ~, typeIdx] = unique(cellTypes);
numTypes = numel(uniqueTypes);

% Initialize count matrix: rows = cell types, columns = direction labels
directionCounts = zeros(numTypes, numDirections);
clusterID = size(cellTypes,1);

% Count occurrences
for i = 1:length(directionLabels)
    row = typeIdx(i);  % cell type index
    col = find(directions == directionLabels(i));  % direction label index
    directionCounts(row, col) = directionCounts(row, col) + 1;
    clusterID(i,1) = col;
end

directionLabels = directionLabels';

eventCounts_clusterID = [eventCounts, clusterID];

      

end
