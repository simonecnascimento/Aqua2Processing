clear all;

% Define the directory containing all experiment files
experimentDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\4._network_propagation.mat\new';

% List all MAT files in the directory
matFiles = dir(fullfile(experimentDir, '*.mat'));

fileNames = {matFiles.name};

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
sortedFileNames = fileNames(idx);

%%
% Initialize an empty table to accumulate results
allNetworkCells = table();
validCells = table();

% Loop through each MAT file
for file = 2:length(sortedFileNames)
    % Load the current experiment's data
    fileName = sortedFileNames{file};
    experimentData = load(fullfile(experimentDir, fileName));

    allCellsNumbersNetwork = [experimentData.networkData.cellNumber_event, experimentData.networkData.cellNumber_network];
    allCellsNumbersNetwork_flattened = cellfun(@(a,b) [a,b], allCellsNumbersNetwork(:,1), allCellsNumbersNetwork(:,2), 'UniformOutput', false);
    allCellsNumbersNetwork_flattened(experimentData.data_analysis.cols_to_delete) = [];

    allCellsTypesNetwork = [experimentData.networkData.cellType_event, experimentData.networkData.cellType_network];
    allCellsTypesNetwork_flattened = cellfun(@(a,b) [a,b], allCellsTypesNetwork(:,1), allCellsTypesNetwork(:,2), 'UniformOutput', false);
    allCellsTypesNetwork_flattened(experimentData.data_analysis.cols_to_delete) = [];

    allCellsNetwork = [allCellsNumbersNetwork_flattened, allCellsTypesNetwork_flattened];
    allCellsNetwork = cell2table(allCellsNetwork, 'VariableNames', {'CellNumbers_all', 'CellTypes_all'});

    allCellsNetwork_cleaned = allCellsNetwork;
    allCellsNetwork_cleaned.InteractionCategory = cell(height(allCellsNetwork_cleaned), 1);
    
    for j = 1:height(allCellsNetwork_cleaned)
        vals1 = allCellsNetwork_cleaned.CellNumbers_all{j};
        vals2 = allCellsNetwork_cleaned.CellTypes_all{j};
    
        % Keep unique cell numbers and align cell types
        [uniqueVals, ia, ~] = unique(vals1, 'stable');
        vals1 = uniqueVals;
        vals2 = vals2(ia);
    
        % Save back cleaned values
        allCellsNetwork_cleaned.CellNumbers_all{j} = vals1;
        allCellsNetwork_cleaned.CellTypes_all{j} = vals2;
    
        % Classify InteractionCategory
        rowVals = unique(vals2);
        if isequal(rowVals, 0)
            allCellsNetwork_cleaned.InteractionCategory{j} = 'P';
        elseif isequal(rowVals, 2)
            allCellsNetwork_cleaned.InteractionCategory{j} = 'NP';
        elseif isequal(rowVals, [0,2])
            allCellsNetwork_cleaned.InteractionCategory{j} = 'mixed';
        else
            allCellsNetwork_cleaned.InteractionCategory{j} = [];
        end
    
        % Remove row if only a single unique value
        if numel(vals1) == 1
            allCellsNetwork_cleaned.CellNumbers_all{j} = [];
            allCellsNeallCellsNetwork_cleanedtwork.CellTypes_all{j}   = [];
            allCellsNetwork_cleaned.InteractionCategory{j} = [];
        end     
    end
        
    % Prepare arrays for expanded table
    interaction_cellIDs = [];
    interaction_categories = [];
    
    % Expand cellNumber_all into rows
    for k = 1:height(allCellsNetwork_cleaned)
        if ~isempty(allCellsNetwork_cleaned.CellNumbers_all{k}) && ~isempty(allCellsNetwork_cleaned.InteractionCategory{k})
            ids = allCellsNetwork_cleaned.CellNumbers_all{k};
            catStr = allCellsNetwork_cleaned.InteractionCategory{k};
            cats = repmat({catStr}, numel(ids), 1);
            
            interaction_cellIDs = [interaction_cellIDs; ids(:)];
            interaction_categories = [interaction_categories; cats];
        end
    end
    
    % Create table of one row per cell
    expandedTable = table(interaction_cellIDs, interaction_categories, ...
        'VariableNames', {'CellID', 'InteractionCategory'});
    
    % Count occurrences by CellID + InteractionCategory
    network_cells = groupsummary(expandedTable, {'CellID', 'InteractionCategory'}, @numel);
    network_cells.Properties.VariableNames{'GroupCount'} = 'Occurrences';

    % Create cellType column (0 = perivascular, 2 = non-perivascular)
    network_cells.cellType = 2 * ones(height(network_cells), 1); % default 2
    isPeri = ismember(network_cells.CellID, experimentData.data_analysis.perivascularCells);
    network_cells.cellType(isPeri) = 0;
    
    % Reorder columns: CellID, cellType, InteractionCategory, Occurrences
    network_cells = network_cells(:, {'CellID', 'cellType', 'InteractionCategory', 'Occurrences'});
    
    % Add a new column to identify the experiment
    network_cells.ExperimentID = repmat({fileName}, height(network_cells), 1);
    
    numberValidCells = numel(experimentData.data_analysis.validCells);
    % Create table
    validCellsTable = table({fileName}, numberValidCells, ...
                            'VariableNames', {'ExperimentID', 'NumValidCells'});

    % Append the current experiment's data to the accumulated table
    allNetworkCells = [allNetworkCells; network_cells];
    validCells = [validCells; validCellsTable];
end

% Save the combined table to a MAT file
save('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\all_network_cells.mat', 'allNetworkCells');

%%
% Example: myTable with columns: cellID, ..., fileName

% Extract FOV names
FOV_names = allNetworkCells.ExperimentID;  % adjust column name if different

% Find unique FOVs
uniqueFOVs = unique(FOV_names);

% Initialize results
cellCounts = zeros(length(uniqueFOVs),1);

% Loop through each FOV and count unique cells
for file = 1:length(uniqueFOVs)
    idx = strcmp(FOV_names, uniqueFOVs{file});
    % Count unique cell IDs in this FOV
    cellCounts(file) = numel(unique(allNetworkCells.CellID(idx)));
end

% Combine into a table for display
FOV_cellCount_table = table(uniqueFOVs, cellCounts, ...
    'VariableNames', {'FOV', 'NumCells'});

disp(FOV_cellCount_table)
sum(FOV_cellCount_table.NumCells)




%%
% Ensure there is an 'InteractionCategory' column for empty rows
% If there are any cells not in the network_cells table, assign 'none'
allNetworkCells.InteractionCategory(cellfun(@isempty, allNetworkCells.InteractionCategory)) = {'none'};

% Define unique types and categories
cellTypes = [0, 2];
categories = {'P','NP','mixed','none'};

% Initialize the 2x4 table
contingencyCounts = zeros(length(cellTypes), length(categories));

% Loop through each cell type and category to count occurrences
for file = 1:length(cellTypes)
    for j = 1:length(categories)
        contingencyCounts(file,j) = sum(allNetworkCells.cellType == cellTypes(file) & ...
                                      strcmp(allNetworkCells.InteractionCategory, categories{j}));
    end
end

% Convert to table with labels
contingencyTable = array2table(contingencyCounts, ...
    'VariableNames', categories, ...
    'RowNames', {'Perivascular_0','NonPerivascular_2'});

disp(contingencyTable);
