% Extract cell types (column 13)
cellTypes = combinedTable_complete{:, 13};

% Unique cell types
[uniqueTypes, ~, typeIdx] = unique(cellTypes);

% Group labels
groupLabels = ["↑↑", "↑→", "↑↓", ...
               "→↑", "→→", "→↓", ...
               "↓↑", "↓→", "↓↓"];

% Initialize result matrix
groupMatrix = zeros(numel(uniqueTypes), 9);

% Fill matrix
for i = 1:numel(uniqueTypes)
    % Group indices for this cell type
    groupForType = groupIndex25(typeIdx == i);
    
    % Count occurrences of each group (1–9)
    counts = histcounts(groupForType, 1:10);
    
    % Store in matrix
    groupMatrix(i, :) = counts;
end

% Create table for display
resultTable = array2table(groupMatrix, 'VariableNames', groupLabels);
resultTable.CellType = uniqueTypes;

% Move CellType column to front
resultTable = movevars(resultTable, 'CellType', 'Before', 1);

% Display
disp(resultTable);


% All possible direction labels
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

% Convert to table for readability
directionTable = array2table(directionCounts, ...
    'VariableNames', directions, ...
    'RowNames', cellstr(uniqueTypes));
