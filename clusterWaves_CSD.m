function [directions, directionCounts, directionLabels, clusterID, rates_mHz, rates_mHz_clusterID] = clusterWaves_CSD(eventHz_byCell, cellTypes)
%CLASSIFYEVENTRATES Classifies cells into 9 groups based on event rate changes.
%   Input:
%     eventHz_byCell - Nx4 string array where:
%       Column 1: Cell ID
%       Column 2: preCSD event rate (Hz)
%       Column 3: duringCSD event rate (Hz)
%       Column 4: postCSD event rate (Hz)
%
%   Output:
%     groupIndex - Nx1 array of group labels (1 to 9)
%     rates_uHz  - Nx3 matrix of event rates in µHz [pre during post]

    % Convert string array to string (in case it's char or mixed)
    eventStr = string(eventHz_byCell);

    % Extract and convert columns
    preHz = str2double(eventStr(:, 2));
    duringHz = str2double(eventStr(:, 3));
    postHz = str2double(eventStr(:, 4));

    % Convert to mHz
    ratesHz = [preHz, duringHz, postHz];
    rates_mHz = ratesHz * 1e3;

    % PRE-TO-DURING Logical conditions for classification of events 
    increasingEventsChange1 = rates_mHz(:, 2) > 2 * rates_mHz(:, 1);
    decreasingEventsChange1 = rates_mHz(:, 2) < 0.5 * rates_mHz(:, 1);
    noChangeEventsChange1 = ~(increasingEventsChange1 | decreasingEventsChange1);

    % PRE-TO-POST Logical conditions for classification of events 
    increasingEventsChange2 = rates_mHz(:, 3) > 2 * rates_mHz(:, 1);
    decreasingEventsChange2 = rates_mHz(:, 3) < 0.5 * rates_mHz(:, 1);
    noChangeEventsChange2 = ~(increasingEventsChange2 | decreasingEventsChange2);

    labels1 = repmat("no change", size(rates_mHz, 1), 1);
    labels1(increasingEventsChange1) = "increasing";
    labels1(decreasingEventsChange1) = "decreasing";

    labels2 = repmat("no change", size(rates_mHz, 1), 1);
    labels2(increasingEventsChange2) = "increasing";
    labels2(decreasingEventsChange2) = "decreasing";

    % Preallocate
    directionLabels = strings(size(rates_mHz,1), 1)';
    
    % Create 9 clusters based on change1 and change2
    for i = 1:length(rates_mHz)
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
        elseif acute == "no change"
            if chronic == "increasing"
                directionLabels(i) = "x↑";
            elseif chronic == "no change"
                directionLabels(i) = "xx";
            elseif chronic == "decreasing"
                directionLabels(i) = "x↓";
            end
        elseif acute == "decreasing"
            if chronic == "increasing"
                directionLabels(i) = "↓↑";
            elseif chronic == "no change"
                directionLabels(i) = "↓x";
            elseif chronic == "decreasing"
                directionLabels(i) = "↓↓";
            end
        end
    end

    directions = ["↑↑", "↑x", "↑↓", "x↑", "xx", "x↓", "↓↑", "↓x", "↓↓"];
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

    rates_mHz_clusterID = [rates_mHz, clusterID];

end

%     Group	pre→during	during→post
%     1	 ↑	↑
%     2	 ↑	x
%     3	 ↑	↓
%     4	 x	↑
%     5	 x	x
%     6	 x	↓
%     7	 ↓	↑
%     8	 ↓	x
%     9	 ↓	↓


