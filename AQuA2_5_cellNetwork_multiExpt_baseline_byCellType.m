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
for file = 1:length(sortedFileNames)
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
    
    for k = 1:height(allCellsNetwork_cleaned)
        ids   = allCellsNetwork_cleaned.CellNumbers_all{k};
        types = allCellsNetwork_cleaned.CellTypes_all{k};
    
        if isempty(ids), continue; end
    
        for j = 2:numel(ids)
            thisID   = ids(j);
            thisType = types(j);
    
            % Partner types = everyone except this cell
            partnerTypes = types(setdiff(1:numel(types), j));
    
            if isempty(partnerTypes)
                category = 'none';
            elseif all(partnerTypes == 0)
                category = 'P';   % only partners are perivascular
            elseif all(partnerTypes == 2)
                category = 'NP';  % only partners are non-perivascular
            else
                category = 'mixed'; % partners include both
            end
    
            % Append
            interaction_cellIDs   = [interaction_cellIDs; thisID];
            interaction_categories = [interaction_categories; {category}];
        end
    end

    % Create expanded table (one row per cell per event)
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

%% Save the combined table to a MAT file
save('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\all_network_cells.mat', '-v7.3');

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
