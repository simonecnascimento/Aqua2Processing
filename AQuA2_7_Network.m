clear; close all; clc;

% Set your experiment folder
networkExperimentDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\4._network_propagation.mat\new';

% Find all *_network_propagation.mat files
FilesAll = dir(fullfile(networkExperimentDir, '*_network_propagation.mat'));

% Extract file names
fileNames = {FilesAll.name};

% Initialize an array to store the extracted parts for sorting
sortKey = [];

% Loop through all filenames to extract the parts for sorting
for file = 1:length(fileNames)
    filename = fileNames{file};
    
    % Extract the number after "Pf4Ai162-" (e.g., '2' from 'Pf4Ai162-2')
    numberAfterPrefix = sscanf(filename, 'Pf4Ai162-%d', 1);
    
    % Extract the date (e.g., '221130' from 'Pf4Ai162-2_221130_FOV6')
    dateStr = regexp(filename, '\d{6}', 'match', 'once');
    
    % Extract the FOV number (e.g., 'FOV6' from 'Pf4Ai162-2_221130_FOV6')
    fovNumber = sscanf(filename, 'Pf4Ai162-%*d_%*d_FOV%d', 1);
    
    % Store the extracted values in a matrix for sorting
    sortKey = [sortKey; numberAfterPrefix, str2double(dateStr), fovNumber];
end

% Sort by the three columns: numberAfterPrefix, date, fovNumber
[~, idx] = sortrows(sortKey);
sortedFiles = FilesAll(idx);
sortedFileNames = fileNames(idx);

% Preallocate results
experimentNames = sortedFileNames;
numCellsAll = zeros(length(sortedFiles), 1);
numCellsAll_initial = zeros(length(sortedFiles), 1);
eventsToDelete_all = cell(length(sortedFiles), 1);
missingCellsNetworkAll = cell(length(sortedFiles), 1);
numRemainingCells = zeros(length(sortedFiles), 1);
percentRemaining = zeros(length(sortedFiles), 1);
numMissingCellsAll = cell(length(sortedFiles), 1);

%% Loop through sorted files
for i = 1:length(sortedFiles)
    filePath = fullfile(sortedFiles(i).folder, sortedFiles(i).name);
    fprintf('Loading file %d/%d: %s\n', i, length(sortedFiles), sortedFiles(i).name);
    
    % Load
    S = load(filePath);

    numCellsAll_initial(i) = S.numCells; % updated total cells for this FOV
    
    % Clean per event
    cell_clean = cell(S.numEvents, 1);
    for j = 1:S.numEvents
        v1 = S.networkData.cellNumber_event{j};         % cells to remove

        v5 = S.networkData.cellNumber_all_unique{j};    % current list
        cell_clean{j} = setdiff(v5(:)', v1);            % cleaned
    end

    eventsToDelete = S.data_analysis.cols_to_delete;
    eventsToDelete_all{i} = eventsToDelete;
    allCells_initial = 1:S.numCells;
    
    cellsToExclude = [];
    for cellIdx = allCells_initial
        eventsInCell = S.data_CFU.cfuInfo1{cellIdx, 2}; % events this cell belongs to
        if all(ismember(eventsInCell, S.data_analysis.cols_to_delete))
            cellsToExclude(end+1) = cellIdx;
        end
    end
  
    % Remove cellToExclude from cell_clean, so is not part of the simultaneous cells
    for j = 1:S.numEvents
        cell_clean{j} = setdiff(cell_clean{j}, cellsToExclude);
    end

    cell_clean(eventsToDelete) = [];

    allCells = setdiff(allCells_initial, cellsToExclude);
    numCellsAll(i) = numel(allCells); % updated total cells for this FOV

    % Flatten all cells inside cell_clean into a single array
    allCellsInCellClean = unique([cell_clean{:}]);
    
    % Check which cells are missing (cells present in allCells but not in cell_clean)
    missingCells = setdiff(allCells, allCellsInCellClean);
    numMissingCells = numel(missingCells);

    % Store missing cells and counts
    missingCellsNetworkAll{i} = missingCells;
    numMissingCellsAll{i} = numMissingCells; % store as numeric for easy plotting


    % Compute remaining and percentage
    numRemainingCells(i) = numel(allCells) - numMissingCells;
    percentRemaining(i) = (numRemainingCells(i) / numel(allCells)) * 100;
    
    % Display status
    if numMissingCells == 0
        disp('✅ All cells have connectivity.');
    else
        fprintf('❌ %d cells are missing connectivity (%.2f%% remaining).\n', ...
            numMissingCells, percentRemaining(i));
        disp('Missing cells:');
        disp(missingCells);
    end

end

% Create and display summary table
T = table(experimentNames', ...
    numCellsAll_initial, ...
    eventsToDelete_all, ...
    cellsToExclude,...
    numCellsAll, ...
    missingCellsNetworkAll, ...
    numMissingCellsAll, ...
    numRemainingCells, ...
    percentRemaining, ...
    'VariableNames', {'ExperimentName', 'TotalCellsInitial', 'EventsToDelete', 'CellsToExclude','TotalCells', 'MissingCellsNetwork', 'NumMissingCells', 'RemainingCells', 'PercentRemaining'});
