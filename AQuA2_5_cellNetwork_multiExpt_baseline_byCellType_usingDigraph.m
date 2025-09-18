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
percentageTable = table();
activeCellsTable = table();
cellCategoryTable = table();

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

    % Make adjacency binary and symmetric
    B = (A > 0);
    B = B | B';

    % Initialize counters
    PP_edges = 0; NP_NP_edges = 0; P_NP_edges = 0;
    PP_cells = []; NP_NP_cells = []; P_NP_cells = [];

    for r = 1:n
        for c = r+1:n
            if B(r,c) && ismember(r,activeCells) && ismember(c,activeCells)
                if strcmp(cellLabels{r},'P') && strcmp(cellLabels{c},'P')
                    PP_edges = PP_edges + 1;
                    PP_cells = [PP_cells, r, c];
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
    PP_cells = unique(PP_cells);
    NP_NP_cells = unique(NP_NP_cells);
    P_NP_cells = unique(P_NP_cells);

    % Total active cells in experiment
    totalActiveCells = numel(activeCells);
    activeCellsTable = [activeCellsTable; table(totalActiveCells)];

    % ---- Classify cells into exclusive categories ----
    onlyP      = setdiff(PP_cells, union(NP_NP_cells, P_NP_cells));
    onlyNP     = setdiff(NP_NP_cells, union(PP_cells, P_NP_cells));
    onlyP_NP   = setdiff(P_NP_cells, union(PP_cells, NP_NP_cells));
    both       = intersect(union(PP_cells, NP_NP_cells), P_NP_cells);
    none       = setdiff(activeCells, union(union(PP_cells, NP_NP_cells), P_NP_cells));

    % Percentages (relative to total active cells)
    pct_onlyP    = numel(onlyP)    / totalActiveCells * 100;
    pct_onlyNP   = numel(onlyNP)   / totalActiveCells * 100;
    pct_onlyP_NP = numel(onlyP_NP) / totalActiveCells * 100;
    pct_both     = numel(both)     / totalActiveCells * 100;
    pct_none     = numel(none)     / totalActiveCells * 100;

    % ---- Store edges counts ----
    interactionsTable = [interactionsTable; 
        table({fileName}, PP_edges, NP_NP_edges, P_NP_edges, ...
        'VariableNames', {'ExperimentID','PP','NP_NP','P_NP'})];

    % ---- Store counts in cellsTable ----
    cellsTable = [cellsTable;
        table({fileName}, numel(PP_cells), numel(NP_NP_cells), numel(P_NP_cells), ...
        numel(both), numel(none), ...
        'VariableNames', {'ExperimentID','PP','NP_NP','P_NP','both','none'})];

    % ---- Store percentages in percentageTable ----
    percentageTable = [percentageTable;
        table({fileName}, pct_onlyP, pct_onlyNP, pct_onlyP_NP, pct_both, pct_none, ...
        'VariableNames', {'ExperimentID','onlyP','onlyNP','onlyP_NP','both','none'})];

end

sum(cellsTable.PP)
sum(cellsTable.NP_NP)
sum(cellsTable.P_NP)
sum(activeCellsTable.totalActiveCells)
