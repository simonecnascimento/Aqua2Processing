clear; close all; clc;

% Set your experiment folder
networkExperimentDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\4._network_propagation.mat';

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
missingCellsAll = cell(length(sortedFiles), 1);
numRemainingCells = zeros(length(sortedFiles), 1);
percentRemaining = zeros(length(sortedFiles), 1);
numMissingCellsAll = cell(length(sortedFiles), 1);

% Loop through sorted files
for i = 1:length(sortedFiles)
    filePath = fullfile(sortedFiles(i).folder, sortedFiles(i).name);
    fprintf('Loading file %d/%d: %s\n', i, length(sortedFiles), sortedFiles(i).name);
    
    % Load
    S = load(filePath);
    
    % Total cells
    allCells = 1:S.numCells;
    numCellsAll(i) = S.numCells;
    
    % Clean per event
    cell_clean = cell(S.numEvents, 1);
    for j = 1:S.numEvents
        v1 = S.networkData.cellNumber_event{j};         % cells to remove

        v5 = S.networkData.cellNumber_all_unique{j};    % current list
        cell_clean{j} = setdiff(v5(:)', v1);            % cleaned
    end
    
    % Unique synchronized cells after removal
    allSyncCells = unique([cell_clean{:}]);
    
    cellNumber_allToDelete = [];
    cellNumber_allToDelete_unique = [];
    for events = 1:size(S.data_analysis.eventsToDelete,2)
        event = S.data_analysis.eventsToDelete(events);
        cellNumber = find(cellfun(@(c) any(c == event), S.data_CFU.cfuInfo1(:, 2)));
        cellNumber_allToDelete = [cellNumber_allToDelete; cellNumber];
        cellNumber_allToDelete_unique = unique(cellNumber_allToDelete);
    end

    allSyncCellsFinal = setdiff(allSyncCells, cellNumber_allToDelete_unique);

   % Missing cells
    missingCells = setdiff(allCells, allSyncCellsFinal);
    missingCellsAll{i} = missingCells;
    numMissingCells = size(missingCells,2);
    numMissingCellsAll{i} = numMissingCells;
    
    % Remaining and percentage
    numRemainingCells(i) = S.numCells - numel(missingCells);
    percentRemaining(i) = (numRemainingCells(i) / S.numCells) * 100;
end

% Create and display summary table
T = table(experimentNames', ...
    numCellsAll, ...
    missingCellsAll, ...
    numMissingCellsAll, ...
    numRemainingCells, ...
    percentRemaining, ...
    'VariableNames', {'ExperimentName', 'TotalCells', 'MissingCells', 'NumMissingCells', 'RemainingCells', 'PercentRemaining'});
