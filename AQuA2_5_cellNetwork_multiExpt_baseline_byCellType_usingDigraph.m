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

% Preallocate tables
interactionsTable = table();
cellsTable = table();
activeCellsTable = table();
percentageTable = table();
percentageTable_P = table();     % percentages relative to P cells
percentageTable_NP = table();    % percentages relative to NP cells


for file = 1:length(sortedFileNames)
    % Load the current experiment's data
    fileName = sortedFileNames{file};
    experimentData = load(fullfile(experimentDir, fileName));

    % Adjacency matrix
    A = experimentData.adjMatrix;  
    n = size(A,1);

    % Cell labels (default NP)
    cellLabels = repmat({'NP'}, n, 1);
    cellLabels(experimentData.data_analysis.perivascularCells) = {'P'};

    % --- Only consider cells with events ---
    cols_to_delete = experimentData.data_analysis.cols_to_delete; % vector of numbers
    col2 = experimentData.event_cell_Table(:,2);                  % get column 2
    
    for i = 1:numel(col2)
        if ~isempty(col2{i})
            % Remove numbers in cols_to_delete from this cell
            col2{i} = col2{i}(~ismember(col2{i}, cols_to_delete));
        end
    end
    
    % Put the cleaned column back into event_cell_Table
    experimentData.event_cell_Table(:,2) = col2;
    activeCells = find(~cellfun(@isempty, experimentData.event_cell_Table(:,2)));
    % Total active cells in experiment
    totalActiveCells = numel(activeCells);
    activeCellsTable = [activeCellsTable; table(totalActiveCells)];

    % Types of cells
    P_cells  = intersect(activeCells, experimentData.data_analysis.perivascularCells);
    NP_cells = setdiff(activeCells, experimentData.data_analysis.perivascularCells);

    % Make adjacency binary and symmetric
    B = (A > 0);
    B = B | B';

    % Initialize counters
    PP_edges = 0; NP_NP_edges = 0; P_NP_edges = 0;
    P_P_cells = []; NP_NP_cells = []; P_NP_cells = [];

    for r = 1:n
        for c = r+1:n
            if B(r,c) && ismember(r,activeCells) && ismember(c,activeCells)
                if strcmp(cellLabels{r},'P') && strcmp(cellLabels{c},'P')
                    PP_edges = PP_edges + 1;
                    P_P_cells = [P_P_cells, r, c];
                elseif strcmp(cellLabels{r},'NP') && strcmp(cellLabels{c},'NP')
                    NP_NP_edges = NP_NP_edges + 1;
                    NP_NP_cells = [NP_NP_cells, r, c];
                else
                    P_NP_edges = P_NP_edges + 1;
                    P_NP_cells = [P_NP_cells, r, c];
                end
            end
        end
    end

    % Unique cells for each interaction type
    P_P_cells = unique(P_P_cells);
    NP_NP_cells = unique(NP_NP_cells);
    P_NP_cells = unique(P_NP_cells);

    % ---- Classify cells into exclusive categories ----
    onlyP      = setdiff(P_P_cells, union(NP_NP_cells, P_NP_cells));
    onlyNP     = setdiff(NP_NP_cells, union(P_P_cells, P_NP_cells));
    onlyP_NP   = setdiff(P_NP_cells, union(P_P_cells, NP_NP_cells));
    %both       = intersect(union(P_P_cells, NP_NP_cells), P_NP_cells);
    %none       = setdiff(activeCells, union(union(P_P_cells, NP_NP_cells), P_NP_cells));

    % "Both" = cells that appear in (P_P or NP_NP) AND in P_NP
    both       = intersect(union(P_P_cells, NP_NP_cells), P_NP_cells);
    
    % Remove those from the others to make categories strictly exclusive
    onlyP      = setdiff(onlyP, both);
    onlyNP     = setdiff(onlyNP, both);
    onlyP_NP   = setdiff(onlyP_NP, both);
    
    allCategorized = union(union(union(onlyP, onlyNP), onlyP_NP), both);

    % Active cells not in any category = none
    none = setdiff(activeCells, allCategorized);

    % ---- Store edges counts ----
    interactionsTable = [interactionsTable; 
        table({fileName}, PP_edges, NP_NP_edges, P_NP_edges, ...
        'VariableNames', {'ExperimentID','PP','NP_NP','P_NP'})];
   
    % Initialize percentages as NaN
    pct_onlyP    = NaN;
    pct_onlyNP   = NaN;
    pct_onlyP_NP = NaN;
    pct_both     = NaN;
    pct_none     = NaN;

    % Compute percentages only if there are active cells
    if totalActiveCells > 0
        if ~isempty(P_cells)
            if numel(P_cells) > 1
                pct_onlyP = numel(onlyP) / totalActiveCells * 100;
            end
        end
        if ~isempty(NP_cells)
            if numel(NP_cells) > 1
                pct_onlyNP = numel(onlyNP) / totalActiveCells * 100;
            end
        end
        if ~isempty(P_cells) && ~isempty(NP_cells)
            pct_onlyP_NP = numel(onlyP_NP) / totalActiveCells * 100;
        end
    
        % both and none can always be computed (they don’t depend on only P or NP existing alone)
        if ~isempty(P_cells) && ~isempty(NP_cells)
           pct_both = numel(both) / totalActiveCells * 100;
        end
        pct_none = numel(none) / totalActiveCells * 100;
    end

    % ---- Store counts in cellsTable ----
    cellsTable = [cellsTable;
        table({fileName}, numel(onlyP), numel(onlyNP), numel(onlyP_NP), ...
        numel(both), numel(none), numel(P_cells), numel(NP_cells), ...
        'VariableNames', {'ExperimentID','PP','NP_NP','P_NP','both','none', 'numP','numNP'})];

    % ---- Store percentages in percentageTable ----
    percentageTable = [percentageTable;
    table({fileName}, pct_onlyP, pct_onlyNP, pct_onlyP_NP, pct_both, pct_none, ...
    numel(P_cells), numel(NP_cells), ...
    'VariableNames', {'ExperimentID','onlyP','onlyNP','onlyP_NP','both','none','numP','numNP'})];

    % Check sum
    sumPct = pct_onlyP + pct_onlyNP + pct_onlyP_NP + pct_both + pct_none;
    if abs(sumPct-100) > 1e-10
        warning('Percentages do not sum to 100 for file %s (sum = %.2f%%)', fileName, sumPct);
    end

%     % ---- Percentages relative to P and NP cells ----
%     % Only P percentages
%     pct_onlyP_relP = numel(onlyP)/max(numel(P_cells),1); % avoid div0
%     pct_both_relP  = numel(both)/max(numel(P_cells),1);
%     pct_none_relP  = numel(setdiff(P_cells, union(onlyP,both)))/max(numel(P_cells),1);
% 
%     percentageTable_P = [percentageTable_P;
%         table({fileName}, pct_onlyP_relP*100, pct_both_relP*100, pct_none_relP*100, ...
%         'VariableNames', {'ExperimentID','onlyP','both','none'})];
% 
%     % Only NP percentages
%     pct_onlyNP_relNP = numel(onlyNP)/max(numel(NP_cells),1);
%     pct_both_relNP   = numel(both)/max(numel(NP_cells),1);
%     pct_none_relNP   = numel(setdiff(NP_cells, union(onlyNP,both)))/max(numel(NP_cells),1);
% 
%     percentageTable_NP = [percentageTable_NP;
%         table({fileName}, pct_onlyNP_relNP*100, pct_both_relNP*100, pct_none_relNP*100, ...
%         'VariableNames', {'ExperimentID','onlyNP','both','none'})];
end

%%

filtered_percentageTable = percentageTable([],:); % keep table structure but empty
for fov = 1:height(percentageTable)
    if percentageTable.numP(fov) > 1 && percentageTable.numNP(fov) > 1
        includedFOV_percentage = percentageTable(fov, :);
        filtered_percentageTable = [filtered_percentageTable; includedFOV_percentage];
    end
end

filtered_cellTable = cellsTable([],:); % keep table structure but empty
for fov = 1:height(cellsTable)
    if cellsTable.numP(fov) > 1 && cellsTable.numNP(fov) > 1
        includedFOV = cellsTable(fov, :);
        filtered_cellTable = [filtered_cellTable; includedFOV];
    end
end


filtered_cellTable_NP = cellsTable([],:); % keep table structure but empty
for fov = 1:height(cellsTable)
    if cellsTable.numP(fov) == 0
        includedFOV_NP = cellsTable(fov, :);
        filtered_cellTable_NP = [filtered_cellTable_NP; includedFOV_NP];
    end
end

filtered_percentageTable_NP = percentageTable([],:); % keep table structure but empty
for fov = 1:height(percentageTable)
    if percentageTable.numP(fov) == 0
        includedFOV_percentage_NP = percentageTable(fov, :);
        filtered_percentageTable_NP = [filtered_percentageTable_NP; includedFOV_percentage_NP];
    end
end


sum(percentageTable.numP)
sum(percentageTable.numNP)
sum(cellsTable.PP)
sum(cellsTable.NP_NP)
sum(cellsTable.P_NP)
sum(activeCellsTable.totalActiveCells)
