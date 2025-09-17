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
validCellsListAll = table();

% Loop through each MAT file
for file = 3:length(sortedFileNames)
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
    
    % Append the current experiment's data to the accumulated table
    allNetworkCells = [allNetworkCells; network_cells];

    % List of valid cells
    % Valid cells - counting also multinucleated here
    numberValidCells_withMN = numel(experimentData.data_analysis.resultsFinal.("Cell ID"));
    validCellsTable = table({fileName}, numberValidCells_withMN, 'VariableNames', {'ExperimentID', 'NumValidCells'});
    validCells = [validCells; validCellsTable];

    validCellsList = experimentData.data_analysis.resultsFinal.("Cell ID");
    validCellsListTable = table(repmat({fileName}, numberValidCells_withMN,1 ), validCellsList, 'VariableNames', {'ExperimentID', 'ValidCells'});
    validCellsListAll = [validCellsListAll; validCellsListTable];
end

sum(validCells.NumValidCells)
% Save the combined table to a MAT file
save('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\all_network_cells.mat');

%%
% Example: myTable with columns: cellID, ..., fileName

% Extract FOV names
FOV_names = validCells.ExperimentID;

% Initialize results
cellCounts = zeros(length(uniqueFOVs),1);

% Loop through each FOV and count unique cells
for file = 1:length(FOV_names)
    idx = strcmp(allNetworkCells.ExperimentID, FOV_names{file});
    % Count unique cell IDs in this FOV
    cellCounts(file) = numel(unique(allNetworkCells.CellID(idx)));
end

% Combine into a table for display
FOV_cellCount_table = table(FOV_names, cellCounts, ...
    'VariableNames', {'FOV', 'NumCells'});

sum(FOV_cellCount_table.NumCells)

%%

% Ensure missing InteractionCategory entries are "none"
allNetworkCells.InteractionCategory(cellfun(@isempty, allNetworkCells.InteractionCategory)) = {'none'};

% Unique cell IDs
allNetworkCells.UniqueCellID = strcat(allNetworkCells.ExperimentID, "_", string(allNetworkCells.CellID));

% Define cell types
cellTypes = [0, 2];

% Initialize summary table with extra columns
summaryCounts = table('Size',[length(cellTypes),8], ...
                      'VariableTypes', repmat("double",1,8), ...
                      'VariableNames', {'OnlyP','OnlyNP','OnlyMixed','NP_or_Mixed','P_or_Mixed','P_and_Mixed','NP_and_Mixed','TotalCells'}, ...
                      'RowNames', {'Perivascular_0','NonPerivascular_2'});

for t = 1:length(cellTypes)
    typeFilter = allNetworkCells.cellType == cellTypes(t);
    
    % Cells per category
    P_cells       = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'P')));
    NP_cells      = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'NP')));
    mixed_cells   = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'mixed')));
    
    % Mutually exclusive counts
    only_P     = setdiff(P_cells, union(NP_cells, mixed_cells));
    only_NP    = setdiff(NP_cells, union(P_cells, mixed_cells));
    only_mixed = setdiff(mixed_cells, union(P_cells, NP_cells));
    
    % Combined categories
    NP_or_mixed  = union(NP_cells, mixed_cells);
    P_or_mixed   = union(P_cells, mixed_cells);
    P_and_mixed  = intersect(P_cells, mixed_cells);  
    NP_and_mixed = intersect(NP_cells, mixed_cells);  % <-- NP and mixed
    
    % Total unique cells of this type
    totalCells = numel(unique(allNetworkCells.UniqueCellID(typeFilter)));
    
    % Save in table
    summaryCounts{t,:} = [numel(only_P), numel(only_NP), numel(only_mixed), ...
                          numel(NP_or_mixed), numel(P_or_mixed), numel(P_and_mixed), numel(NP_and_mixed), totalCells];
end

disp(summaryCounts)


%%
% Initialize summary table
summaryCounts = table('Size',[length(cellTypes),6], ...
                      'VariableTypes', repmat("double",1,6), ...
                      'VariableNames', {'OnlyP','OnlyNP','OnlyMixed','NP_or_Mixed','P_or_Mixed','TotalCells'}, ...
                      'RowNames', {'Perivascular_0','NonPerivascular_2'});

for t = 1:length(cellTypes)
    typeFilter = allNetworkCells.cellType == cellTypes(t);
    
    % Cells per category
    P_cells       = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'P')));
    NP_cells      = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'NP')));
    mixed_cells   = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'mixed')));
    
    % Mutually exclusive counts
    only_P     = setdiff(P_cells, union(NP_cells, mixed_cells));
    only_NP    = setdiff(NP_cells, union(P_cells, mixed_cells));
    only_mixed = setdiff(mixed_cells, union(P_cells, NP_cells));
    
    % Combined categories
    NP_or_mixed = union(NP_cells, mixed_cells);
    P_or_mixed  = union(P_cells, mixed_cells);
    
    % Total unique cells of this type
    totalCells = numel(unique(allNetworkCells.UniqueCellID(typeFilter)));
    
    % Save in table
    summaryCounts{t,:} = [numel(only_P), numel(only_NP), numel(only_mixed), ...
                          numel(NP_or_mixed), numel(P_or_mixed), totalCells];
end

disp(summaryCounts)
%%

% Initialize summary table with extra columns
summaryCounts = table('Size',[length(cellTypes),8], ...
                      'VariableTypes', repmat("double",1,8), ...
                      'VariableNames', {'OnlyP','OnlyNP','OnlyMixed','NP_or_Mixed','P_or_Mixed','P_and_Mixed','NP_and_Mixed','TotalCells'}, ...
                      'RowNames', {'Perivascular_0','NonPerivascular_2'});

for t = 1:length(cellTypes)
    typeFilter = allNetworkCells.cellType == cellTypes(t);
    
    % Cells per category
    P_cells       = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'P')));
    NP_cells      = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'NP')));
    mixed_cells   = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'mixed')));
    
    % Mutually exclusive counts
    only_P     = setdiff(P_cells, union(NP_cells, mixed_cells));
    only_NP    = setdiff(NP_cells, union(P_cells, mixed_cells));
    only_mixed = setdiff(mixed_cells, union(P_cells, NP_cells));
    
    % Combined categories
    NP_or_mixed  = union(NP_cells, mixed_cells);
    P_or_mixed   = union(P_cells, mixed_cells);
    P_and_mixed  = intersect(P_cells, mixed_cells);  
    NP_and_mixed = intersect(NP_cells, mixed_cells);  % <-- NP and mixed
    
    % Total unique cells of this type
    totalCells = numel(unique(allNetworkCells.UniqueCellID(typeFilter)));
    
    % Save in table
    summaryCounts{t,:} = [numel(only_P), numel(only_NP), numel(only_mixed), ...
                          numel(NP_or_mixed), numel(P_or_mixed), numel(P_and_mixed), numel(NP_and_mixed), totalCells];
end

disp(summaryCounts)


%%
% Ensure missing InteractionCategory entries are "none"
allNetworkCells.InteractionCategory(cellfun(@isempty, allNetworkCells.InteractionCategory)) = {'none'};

% Unique cell IDs
allNetworkCells.UniqueCellID = strcat(allNetworkCells.ExperimentID, "_", string(allNetworkCells.CellID));

% Define cell types
cellTypes = [0, 2];

% Initialize summary table
summaryCounts = table('Size',[length(cellTypes),5], ...
                      'VariableTypes', repmat("double",1,5), ...
                      'VariableNames', {'OnlyP','OnlyNP','OnlyMixed','NPorMixed','TotalCells'}, ...
                      'RowNames', {'Perivascular_0','NonPerivascular_2'});

for t = 1:length(cellTypes)
    typeFilter = allNetworkCells.cellType == cellTypes(t);
    
    % Cells per category
    P_cells       = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'P')));
    NP_cells      = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'NP')));
    mixed_cells   = unique(allNetworkCells.UniqueCellID(typeFilter & strcmp(allNetworkCells.InteractionCategory,'mixed')));
    
    % Mutually exclusive
    only_P     = setdiff(P_cells, union(NP_cells, mixed_cells));
    only_NP    = setdiff(NP_cells, union(P_cells, mixed_cells));
    only_mixed = setdiff(mixed_cells, union(P_cells, NP_cells));
    NP_or_mixed = union(NP_cells, mixed_cells);
    
    % Total unique cells of this type
    totalCells = numel(unique(allNetworkCells.UniqueCellID(typeFilter)));
    
    % Save in table
    summaryCounts{t,:} = [numel(only_P), numel(only_NP), numel(only_mixed), numel(NP_or_mixed), totalCells];
end

disp(summaryCounts)
%%

allNetworkCells.UniqueCellID = strcat(allNetworkCells.ExperimentID, "_", string(allNetworkCells.CellID));

[uniqueCells, ~, G] = unique(allNetworkCells.UniqueCellID); % uniqueCells = original cell IDs
counts = accumarray(G, 1); % how many times each unique cell appears

cellTypes   = splitapply(@(x) x(1), allNetworkCells.cellType, G); % each cell has only one type
categories  = splitapply(@(x) {unique(x)}, allNetworkCells.InteractionCategory, G);

perCellSummary = table(uniqueCells, cellTypes, categories, ...
    'VariableNames', {'UniqueCellID','cellType','Categories'});

% Ensure missing InteractionCategory entries are "none"
allNetworkCells.InteractionCategory(cellfun(@isempty, allNetworkCells.InteractionCategory)) = {'none'};

% Define types and categories
cellTypes   = [0, 2];
categories  = {'P','NP','mixed','none'};

% Initialize table
contingencyCounts = zeros(length(cellTypes), length(categories));

for t = 1:length(cellTypes)
    for c = 1:length(categories)
        % Count *unique cells* of this type in this category
        uniqueCells = unique(allNetworkCells.UniqueCellID( ...
            allNetworkCells.cellType == cellTypes(t) & ...
            strcmp(allNetworkCells.InteractionCategory, categories{c})));
        contingencyCounts(t,c) = numel(uniqueCells);
    end
end

contingencyTable = array2table(contingencyCounts, ...
    'VariableNames', categories, ...
    'RowNames', {'Perivascular_0','NonPerivascular_2'});


% Cells in NP or mixed
NP_or_mixed = unique(allNetworkCells.UniqueCellID( ...
    ismember(allNetworkCells.InteractionCategory, {'NP','mixed'})));

% Cells only in mixed
only_mixed = setdiff( ...
    unique(allNetworkCells.UniqueCellID(strcmp(allNetworkCells.InteractionCategory,'mixed'))), ...
    unique(allNetworkCells.UniqueCellID(strcmp(allNetworkCells.InteractionCategory,'NP'))));

fprintf('Cells in NP or mixed: %d\n', numel(NP_or_mixed));
fprintf('Cells only in mixed: %d\n', numel(only_mixed));



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

%% --- Additional checks ---

% Total cells in dataset
totalCells = sum(validCells.NumValidCells);

allNetworkCells.UniqueCellID = strcat(allNetworkCells.ExperimentID, "_", string(allNetworkCells.CellID));

% Cells that ever participated in NP or mixed
NP_or_mixed_cells = unique(allNetworkCells.CellID( ...
    ismember(allNetworkCells.InteractionCategory, {'NP','mixed'})));

% Cells that only participated in mixed (not NP)
only_mixed_cells = setdiff(unique(allNetworkCells.CellID(strcmp(allNetworkCells.InteractionCategory,'mixed'))), ...
                           unique(allNetworkCells.CellID(strcmp(allNetworkCells.InteractionCategory,'NP'))));

fprintf('Total unique cells: %d\n', totalCells);
fprintf('Cells in NP or mixed: %d\n', numel(NP_or_mixed_cells));
fprintf('Cells only in mixed: %d\n', numel(only_mixed_cells));
