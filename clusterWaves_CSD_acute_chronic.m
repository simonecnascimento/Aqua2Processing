function [directions, acuteDirectionLabels, chronicDirectionLabels, ...
          acuteClusterID, chronicClusterID, ...
          acuteDirectionCounts, chronicDirectionCounts, ...
          rates_mHz] = clusterWaves_CSD_acute_chronic(eventHz_byCell, cellTypes)
%CLUSTERWAVES_CSD_ACUTE_CHRONIC Classifies cells into ↑ / ↓ / NA for both acute and chronic changes.
%
% Inputs:
%   eventHz_byCell - Nx4 string array: {cellID, preHz, duringHz, postHz}
%   cellTypes      - Nx1 string or categorical array indicating cell types
%
% Outputs:
%   acuteDirectionLabels   - Nx1 string array ("↑", "↓", or "NA") for pre → during
%   chronicDirectionLabels - Nx1 string array ("↑", "↓", or "NA") for pre → post
%   acuteClusterID         - Nx1 numeric (1=↑, 2=↓, NaN=NA)
%   chronicClusterID       - Nx1 numeric (1=↑, 2=↓, NaN=NA)
%   acuteDirectionCounts   - [numCellTypes x 2] matrix: ↑ and ↓ counts by type
%   chronicDirectionCounts - [numCellTypes x 2] matrix: ↑ and ↓ counts by type
%   rates_mHz              - Nx3 matrix: [pre, during, post] in mHz

    % Convert to string and extract rates
    eventStr = string(eventHz_byCell);
    preHz    = str2double(eventStr(:, 2));
    duringHz = str2double(eventStr(:, 3));
    postHz   = str2double(eventStr(:, 4));

    % Convert to mHz
    rates_mHz = [preHz, duringHz, postHz] * 1e3;

    % Define thresholds
    acuteUp   = rates_mHz(:,2) > 2 * rates_mHz(:,1);
    acuteDown = rates_mHz(:,2) < 0.5 * rates_mHz(:,1);
    chronicUp   = rates_mHz(:,3) > 2 * rates_mHz(:,1);
    chronicDown = rates_mHz(:,3) < 0.5 * rates_mHz(:,1);

    % Initialize
    N = size(rates_mHz,1);
    acuteDirectionLabels   = repmat("NA", N, 1);
    chronicDirectionLabels = repmat("NA", N, 1);
    acuteClusterID         = nan(N, 1);
    chronicClusterID       = nan(N, 1);

    % Label directions
    acuteDirectionLabels(acuteUp)   = "↑";
    acuteDirectionLabels(acuteDown) = "↓";
    acuteClusterID(acuteUp)   = 1;
    acuteClusterID(acuteDown) = 2;

    chronicDirectionLabels(chronicUp)   = "↑";
    chronicDirectionLabels(chronicDown) = "↓";
    chronicClusterID(chronicUp)   = 1;
    chronicClusterID(chronicDown) = 2;

    % Count directions per cell type
    directions = ["↑", "↓"];
    [uniqueTypes, ~, typeIdx] = unique(cellTypes);
    numTypes = numel(uniqueTypes);
    numDirections = numel(directions);

    acuteDirectionCounts = zeros(numTypes, numDirections);
    chronicDirectionCounts = zeros(numTypes, numDirections);

    for i = 1:N
        row = typeIdx(i);
        % Acute
        aLabel = acuteDirectionLabels(i);
        if aLabel ~= "NA"
            col = find(directions == aLabel);
            acuteDirectionCounts(row, col) = acuteDirectionCounts(row, col) + 1;
        end
        % Chronic
        cLabel = chronicDirectionLabels(i);
        if cLabel ~= "NA"
            col = find(directions == cLabel);
            chronicDirectionCounts(row, col) = chronicDirectionCounts(row, col) + 1;
        end
    end

    directions = ["↑", "↓"];
    cells = ["Perivascular", "Non-perivascular"];

    %Acute
    figure;
    bar(acuteDirectionCounts,'DisplayName','acuteDirectionCounts')
    xticklabels(cells);
    title("Acute directions by cell type")
    legend(directions, Location="best");

    %Chronic
    figure;
    bar(chronicDirectionCounts,'DisplayName','chronicDirectionCounts')
    xticklabels(cells);
    title("Chronic directions by cell type")
    legend(directions, Location="best");

end
