function [directionLabels, directionCounts, clusterID, rates_mHz] = clusterWaves_CSD_simple(eventHz_byCell, cellTypes, mode)
%CLUSTERWAVES_CSD_SIMPLE Clusters cells into 2 groups (↑ or ↓) based on acute or chronic change.
%
% Inputs:
%   eventHz_byCell - Nx4 string array: {cellID, preHz, duringHz, postHz}
%   cellTypes      - Nx1 string or categorical vector for cell types
%   mode           - 'acute' (pre → during) or 'chronic' (pre → post)
%
% Outputs:
%   directionLabels - Nx1 string: "↑", "↓", or "NA"
%   directionCounts - [numCellTypes x 2] count of ↑ and ↓ per type
%   clusterID       - Nx1 group index: 1 (↑), 2 (↓), NaN if no clear change
%   rates_mHz       - [pre, during, post] in mHz

    % Parse input
    eventStr = string(eventHz_byCell);
    preHz    = str2double(eventStr(:, 2));
    duringHz = str2double(eventStr(:, 3));
    postHz   = str2double(eventStr(:, 4));

    % Convert to mHz
    rates_mHz = [preHz, duringHz, postHz] * 1e3;

    % Choose which direction to classify
    switch lower(mode)
        case 'acute'
            compareTo = rates_mHz(:,2); % during
        case 'chronic'
            compareTo = rates_mHz(:,3); % post
        otherwise
            error('Mode must be "acute" or "chronic"');
    end

    % Classify ↑ or ↓ (using >2x or <0.5x thresholds)
    isUp   = compareTo > 2 * rates_mHz(:,1);
    isDown = compareTo < 0.5 * rates_mHz(:,1);

    % Assign labels
    directionLabels = strings(size(rates_mHz,1), 1);
    directionLabels(isUp)   = "↑";
    directionLabels(isDown) = "↓";
    directionLabels(~(isUp | isDown)) = "NA";

    % Define output directions and initialize counts
    directions = ["↑", "↓"];
    numDirections = numel(directions);
    [uniqueTypes, ~, typeIdx] = unique(cellTypes);
    numTypes = numel(uniqueTypes);
    
    directionCounts = zeros(numTypes, numDirections);
    clusterID = nan(size(cellTypes,1), 1);

    % Count per cell type
    for i = 1:length(directionLabels)
        label = directionLabels(i);
        if label ~= "NA"
            row = typeIdx(i);
            col = find(directions == label);
            directionCounts(row, col) = directionCounts(row, col) + 1;
            clusterID(i) = col;
        end
    end
end
