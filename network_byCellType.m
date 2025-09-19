allCellsNumbersNetwork = [networkData.cellNumber_event, networkData.cellNumber_network];
allCellsNumbersNetwork_flattened = cellfun(@(a,b) [a,b], allCellsNumbersNetwork(:,1), allCellsNumbersNetwork(:,2), 'UniformOutput', false);
allCellsNumbersNetwork_flattened(data_analysis.cols_to_delete) = [];

allCellsTypesNetwork = [networkData.cellType_event, networkData.cellType_network];
allCellsTypesNetwork_flattened = cellfun(@(a,b) [a,b], allCellsTypesNetwork(:,1), allCellsTypesNetwork(:,2), 'UniformOutput', false);
allCellsTypesNetwork_flattened(data_analysis.cols_to_delete) = [];

allCellsNetwork = [allCellsNumbersNetwork_flattened, allCellsTypesNetwork_flattened];


% Find rows with more than 1 element
rowsToKeep = cellfun(@numel, allCellsTypesNetwork_flattened) > 1;
% Keep only those rows
allCellsNetwork_flattened_notSingleRow = allCellsTypesNetwork_flattened(rowsToKeep, :);

numOnly0 = 0;
numOnly2 = 0;
numMixed = 0;

for i = 1:numel(allCellsNetwork_flattened_notSingleRow)
    rowVals = unique(allCellsNetwork_flattened_notSingleRow{i});  % get unique values in this row
    
    if isequal(rowVals, 0)  % just 0’s
        numOnly0 = numOnly0 + 1;
    elseif isequal(rowVals, 2)  % just 2’s
        numOnly2 = numOnly2 + 1;
    elseif all(ismember(rowVals, [0 2]))  % contains both 0 and 2
        numMixed = numMixed + 1;
    end
end

fprintf('total events: %d\n', numel(allCellsTypesNetwork_flattened));
fprintf('Rows with only 0: %d\n', numOnly0);
fprintf('Rows with only 2: %d\n', numOnly2);
fprintf('Rows with both 0 and 2: %d\n', numMixed);
numSingle = sum(cellfun(@numel, allCellsTypesNetwork_flattened) == 1);
fprintf('Number of single-element rows: %d\n', numSingle);

%%

cells_sharedEvents = [allCellsTypesNetwork_flattened, networkData.cellNumber_all_unique];
rowsToKeep = cellfun(@numel, cells_sharedEvents(:,2)) > 1;
cells_sharedEvents_notSingleRow = cells_sharedEvents(rowsToKeep, :);

% Extract all numbers from the 2nd column
allNumbers = [cells_sharedEvents_notSingleRow{:,2}];   % concatenates all vectors

% Count unique numbers
uniqueNums = unique(allNumbers);
counts = histc(allNumbers, uniqueNums);

% Display result
T = table(uniqueNums(:), counts(:), 'VariableNames', {'CellID','Count'});
disp(T)

%%
% Extract the two columns into a new table
network_cellData = networkData(:, {'cellNumber_all', 'cellType_all'});

% Initialize InteractionCategory column
network_cellData.InteractionCategory = cell(height(network_cellData),1);

% Apply rules to classify InteractionCategory
for i = 1:height(network_cellData)
    emptyRow = false;

    % Check cellNumber_all for duplicate cells only [1,1] or single number (= no simultaneous events, only its own event)
    vals1 = network_cellData.cellNumber_all{i};
    if isnumeric(vals1) && (~isscalar(vals1) && numel(unique(vals1)) == 1 || isscalar(vals1))
        emptyRow = true;
    end

    [uniqueVals, ia, ~] = unique(vals1, 'stable');  % keep first occurrence
    vals1 = uniqueVals;
    vals2 = vals2(ia);  % keep corresponding cellType values

    if emptyRow
        network_cellData.cellNumber_all{i} = [];
        network_cellData.cellType_all{i}  = [];
        network_cellData.InteractionCategory{i}  = [];
    else
        vals2 = network_cellData.cellType_all{i};
        if ~isempty(vals2)
            rowVals = unique(vals2);
            if isequal(rowVals, 0)
                network_cellData.InteractionCategory{i} = 'P';
            elseif isequal(rowVals, 2)
                network_cellData.InteractionCategory{i} = 'NP';
            elseif all(ismember(rowVals, [0 2]))
                network_cellData.InteractionCategory{i} = 'mixed';
            end
        else
            network_cellData.InteractionCategory{i} = 'empty';
        end
    end
end

% Prepare arrays for expanded table
interaction_cellIDs   = [];
interaction_categories = [];

% Expand cellNumber_all into rows
for i = 1:height(network_cellData)
    if ~isempty(network_cellData.cellNumber_all{i}) && ~isempty(network_cellData.InteractionCategory{i})
        ids = network_cellData.cellNumber_all{i};
        catStr = network_cellData.InteractionCategory{i};
        cats = repmat({catStr}, numel(ids), 1);

        interaction_cellIDs    = [interaction_cellIDs; ids(:)];
        interaction_categories = [interaction_categories; cats];
    end
end

% Create table of one row per cell
expandedTable = table(interaction_cellIDs, interaction_categories, ...
    'VariableNames', {'CellID','InteractionCategory'});

% Count occurrences by CellID + InteractionCategory
network_cells = groupsummary(expandedTable, {'CellID','InteractionCategory'}, @numel);
network_cells.Properties.VariableNames{'GroupCount'} = 'Occurrences';

% Create cellType column (0 = perivascular, 2 = non-perivascular)
network_cells.cellType = 2 * ones(height(network_cells),1); % default 2
isPeri = ismember(network_cells.CellID, data_analysis.perivascularCells);
network_cells.cellType(isPeri) = 0;

% Reorder columns: CellID, cellType, InteractionCategory, Occurrences
network_cells = network_cells(:, {'CellID','cellType','InteractionCategory','Occurrences'});
