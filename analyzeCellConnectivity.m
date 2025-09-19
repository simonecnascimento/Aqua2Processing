function NetworkTable = analyzeCellConnectivity(numCells, numEvents, networkData, data_analysis, data_CFU)
% analyzeCellConnectivity - Analyze cell connectivity in one FOV
%
% Inputs:
%   numCells      - total number of cells in the FOV
%   numEvents     - total number of events in the FOV
%   networkData   - struct with fields:
%                   .cellNumber_event (cell array per event)
%                   .cellNumber_all_unique (cell array per event)
%   data_analysis - struct with fields:
%                   .cols_to_delete (indices of events to delete)
%                   .fileTemp (string, experiment/FOV name)
%   data_CFU      - struct with field:
%                   .cfuInfo1 (cell array, row = cell, col 2 = events this cell belongs to)
%
% Output:
%   NetworkTable  - struct containing analysis results and summary table T

% Store experiment name
ExperimentName = data_analysis.fileTemp;

% Store initial total cells
numCells = numel(data_analysis.validCells);
NetworkTable.numCellsAll_initial = numCells;

% Clean cells per event by removing cells to exclude
NetworkTable.cell_clean = cell(numEvents, 1);
for j = 1:numEvents
    cellsToRemove = networkData.cellNumber_event{j};
    currentCells = networkData.cellNumber_all_unique{j};
    NetworkTable.cell_clean{j} = setdiff(currentCells(:)', cellsToRemove);
end

% Store events to delete
NetworkTable.eventsToDelete = data_analysis.cols_to_delete;

% List all cells initially
NetworkTable.allCells_initial = 1:numCells;

% Identify cells to exclude (cells whose all events are in delete list)
NetworkTable.cellsToExclude = [];
for cellIdx = NetworkTable.allCells_initial
    eventsThisCell = data_CFU.cfuInfo1{cellIdx, 2};
    if all(ismember(eventsThisCell, NetworkTable.eventsToDelete))
        NetworkTable.cellsToExclude(end+1) = cellIdx;
    end
end

% Remove excluded cells from cell_clean
for j = 1:numEvents
    NetworkTable.cell_clean{j} = setdiff(NetworkTable.cell_clean{j}, NetworkTable.cellsToExclude);
end

% Remove deleted events from cell_clean
NetworkTable.cell_clean(NetworkTable.eventsToDelete) = [];

% Update all cells after exclusion
NetworkTable.allCells = setdiff(NetworkTable.allCells_initial, NetworkTable.cellsToExclude);

% Count cells remaining
NetworkTable.numCellsAll = numel(NetworkTable.allCells);

% Flatten all cleaned cells into unique array
NetworkTable.allCellsInCellClean = unique([NetworkTable.cell_clean{:}]);

% % Find missing cells (present in allCells but missing in cleaned cells)
% NetworkTable.missingCells = setdiff(NetworkTable.allCells, NetworkTable.allCellsInCellClean);
% NetworkTable.numMissingCells = numel(NetworkTable.missingCells);

% Find missing cells (present in allCells but missing in cleaned cells)
rawMissing = setdiff(NetworkTable.allCells, NetworkTable.allCellsInCellClean);

% Exclude cells that were merged into others
if isfield(data_analysis, 'cellsToMerge') && ~isempty(data_analysis.cellsToMerge)
    % data_analysis.cellsToMerge is assumed to be N x 2 array: [sourceCell mergedIntoCell]
    mergedCells = unique(data_analysis.cellsToMerge(:));  % flatten all merged cells
    NetworkTable.missingCells = setdiff(rawMissing, mergedCells);  % only truly missing
else
    NetworkTable.missingCells = rawMissing;
end

NetworkTable.numMissingCells = numel(NetworkTable.missingCells);

% Compute remaining cells and percentage
NetworkTable.numRemainingCells = NetworkTable.numCellsAll - NetworkTable.numMissingCells;
NetworkTable.percentRemaining = (NetworkTable.numRemainingCells / NetworkTable.numCellsAll) * 100;

% Display summary
if NetworkTable.numMissingCells == 0
    disp('✅ All cells have connectivity.');
else
    fprintf('❌ %d cells are missing connectivity (%.2f%% remaining).\n', NetworkTable.numMissingCells, NetworkTable.percentRemaining);
    disp('Missing cells:');
    disp(NetworkTable.missingCells);
end

% Create summary table T
NetworkTable.T = table( ...
    {ExperimentName}, ...
    NetworkTable.numCellsAll_initial, ...
    {NetworkTable.eventsToDelete}, ...
    {NetworkTable.cellsToExclude}, ...
    NetworkTable.numCellsAll, ...
    {NetworkTable.missingCells}, ...
    NetworkTable.numMissingCells, ...
    NetworkTable.numRemainingCells, ...
    NetworkTable.percentRemaining, ...
    'VariableNames', {'ExperimentName', 'TotalCellsInitial', 'EventsToDelete', 'CellsToExclude', ...
                      'TotalCells', 'MissingCellsNetwork', 'NumMissingCells', 'RemainingCells', 'PercentRemaining'} ...
);

end
