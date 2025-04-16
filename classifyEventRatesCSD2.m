function [directions, directionCounts, rates_mHz] = classifyEventRatesCSD(eventHz_byCell, cellTypes)
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

    % Convert to µHz
    ratesHz = [preHz, duringHz, postHz];
    rates_mHz = ratesHz * 1e3;

%     % Thresholds for 25% change
%     upThresh = 1.25;
%     downThresh = 0.75;
%     epsilon = 1e-10;  % avoid division by zero
% 
%     % Phase-to-phase ratios
%     preToDuring = rates_mHz(:,2) ./ (rates_mHz(:,1) + epsilon);
%     preToPost = rates_mHz(:,3) ./ (rates_mHz(:,1) + epsilon);
% 
%     % Change classification: 1 = increase, 2 = no change, 3 = decrease
%     change1 = ones(size(preToDuring)) * 2;
%     change1(preToDuring > upThresh) = 1;
%     change1(preToDuring < downThresh) = 3;
% 
%     change2 = ones(size(preToPost)) * 2;
%     change2(preToPost > upThresh) = 1;
%     change2(preToPost < downThresh) = 3;
% 
%     % Group index: 1 to 9
%     groupIndex25 = (change1 - 1) * 3 + change2;

    % PRE-TO-DURING Logical conditions for classification of events 
    risingEventsChange1 = rates_mHz(:, 2) > 2 * rates_mHz(:, 1);
    decreasingEventsChange1 = rates_mHz(:, 2) < 0.5 * rates_mHz(:, 1);
    noChangeEventsChange1 = ~(risingEventsChange1 | decreasingEventsChange1);

     % PRE-TO-POST Logical conditions for classification of events 
    risingEventsChange2 = rates_mHz(:, 3) > 2 * rates_mHz(:, 1);
    decreasingEventsChange2 = rates_mHz(:, 3) < 0.5 * rates_mHz(:, 1);
    noChangeEventsChange2 = ~(risingEventsChange2 | decreasingEventsChange2);

    labels1 = repmat("no change", size(rates_mHz, 1), 1);
    labels1(risingEventsChange1) = "rising";
    labels1(decreasingEventsChange1) = "decreasing";

    labels2 = repmat("no change", size(rates_mHz, 1), 1);
    labels2(risingEventsChange2) = "rising";
    labels2(decreasingEventsChange2) = "decreasing";

    % Preallocate
    directionLabels = strings(size(rates_mHz,1), 1)';
    
    % Create 9 clusters based on change1 and change2
    for i = 1:length(rates_mHz)
        from = labels1(i);  % pre → during
        to = labels2(i);    % pre → post
    
        if from == "rising"
            if to == "rising"
                directionLabels(i) = "↑↑";
            elseif to == "no change"
                directionLabels(i) = "↑→";
            elseif to == "decreasing"
                directionLabels(i) = "↑↓";
            end
        elseif from == "no change"
            if to == "rising"
                directionLabels(i) = "→↑";
            elseif to == "no change"
                directionLabels(i) = "→→";
            elseif to == "decreasing"
                directionLabels(i) = "→↓";
            end
        elseif from == "decreasing"
            if to == "rising"
                directionLabels(i) = "↓↓";
            elseif to == "no change"
                directionLabels(i) = "↓→";
            elseif to == "decreasing"
                directionLabels(i) = "↓↓";
            end
        end
    end

    directions = ["↑↑", "↑→", "↑↓", "→↑", "→→", "→↓", "↓↓", "↓→", "↓↑"];
    numDirections = numel(directions);

    % Unique cell types
    [uniqueTypes, ~, typeIdx] = unique(cellTypes);
    numTypes = numel(uniqueTypes);
    
    % Initialize count matrix: rows = cell types, columns = direction labels
    directionCounts = zeros(numTypes, numDirections);

    % Count occurrences
    for i = 1:length(directionLabels)
        row = typeIdx(i);  % cell type index
        col = find(directions == directionLabels(i));  % direction label index
        directionCounts(row, col) = directionCounts(row, col) + 1;
    end

end

%     Group	pre→during	during→post
%     1	 ↑	↑
%     2	 ↑	→
%     3	 ↑	↓
%     4	 →	↑
%     5	 →	→
%     6	 →	↓
%     7	 ↓	↑
%     8	 ↓	→
%     9	 ↓	↓


