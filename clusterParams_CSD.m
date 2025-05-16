function [directions, directionCounts, directionLabels, clusterID, phaseParams_dFF] = clusterParams_CSD(paramTables_allPhases, cellTypes)
%CLASSIFYEVENTRATES Classifies cells into 9 groups based on event rate changes.
%   Input:
%     eventHz_byCell - Nx4 string array where:
%       Column 1: Cell ID
%       Column 2: preCSD event rate (Hz)
%       Column 3: duringCSD event rate (Hz)
%       Column 4: postCSD event rate (Hz)
%
%   Output:
%     groupIndex - Nx1 array of group labels (1 to 9)
%     rates_uHz  - Nx3 matrix of event rates in µHz [pre during post]

fields = fieldnames(paramTables_allPhases);  

for x = 1:numel(fieldnames(paramTables_allPhases))
    fieldName = fields{x}; 

    % Extract and convert columns
    param_preCSD = paramTables_allPhases.(fieldName)(:, 2);
    param_duringCSD = paramTables_allPhases.(fieldName)(:, 3);
    param_postCSD = paramTables_allPhases.(fieldName)(:, 4);

    phaseParams = [param_preCSD, param_duringCSD, param_postCSD];

    % Replace empty cells with 0
    for i = 1:width(phaseParams)
        if iscell(phaseParams{:, i})
            emptyRows = cellfun(@(x) isempty(x) || isequal(size(x), [0 0]), phaseParams{:, i});
            phaseParams{emptyRows, i} = {0};
        end
    end

    if x == 4
        phaseParams_dFF = phaseParams;
    end

    % Convert each column to double, handling single or numeric values
    numRows = size(phaseParams, 1);
    numCols = size(phaseParams, 2);
    phaseParamsDouble = zeros(numRows, numCols);
    
    for col = 1:numCols
        colData = phaseParams{:, col}; % get the column data as a cell or numeric array
        
        if iscell(colData)
            % If cell, convert each cell element to double
            for row = 1:numRows
                val = colData{row};
                if isempty(val)
                    phaseParamsDouble(row, col) = 0;  % replace empty with 0
                else
                    phaseParamsDouble(row, col) = double(val);
                end
            end
        else
            % If numeric array (single or double), convert whole column to double
            phaseParamsDouble(:, col) = double(colData);
        end
    end

    % PRE-TO-DURING Logical conditions for classification of events 
    increasingParamChange1 = phaseParamsDouble(:, 2) > 2 * phaseParamsDouble(:, 1);
    decreasingParamChange1 = phaseParamsDouble(:, 2) < 0.5 * phaseParamsDouble(:, 1);
    noChangeParamChange1 = ~(increasingParamChange1 | decreasingParamChange1);

    % PRE-TO-POST Logical conditions for classification of events 
    increasingParamChange2 = phaseParamsDouble(:, 3) > 2 * phaseParamsDouble(:, 1);
    decreasingParamChange2 = phaseParamsDouble(:, 3) < 0.5 * phaseParamsDouble(:, 1);
    noChangeParamChange2 = ~(increasingParamChange2 | decreasingParamChange2);

    labels1 = repmat("no change", size(phaseParams, 1), 1);
    labels1(increasingParamChange1) = "increasing";
    labels1(decreasingParamChange1) = "decreasing";

    labels2 = repmat("no change", size(phaseParams, 1), 1);
    labels2(increasingParamChange2) = "increasing";
    labels2(decreasingParamChange2) = "decreasing";

    % Preallocate
    directionLabels = strings(size(phaseParamsDouble,1), 1)';
    
    % Create 9 clusters based on change1 and change2
    for i = 1:length(phaseParamsDouble)
        acute = labels1(i);  % pre → during
        chronic = labels2(i);    % pre → post
    
        if acute == "increasing"
            if chronic == "increasing"
                directionLabels(i) = "↑↑";
            elseif chronic == "no change"
                directionLabels(i) = "↑x";
            elseif chronic == "decreasing"
                directionLabels(i) = "↑↓";
            end
        elseif acute == "no change"
            if chronic == "increasing"
                directionLabels(i) = "x↑";
            elseif chronic == "no change"
                directionLabels(i) = "xx";
            elseif chronic == "decreasing"
                directionLabels(i) = "x↓";
            end
        elseif acute == "decreasing"
            if chronic == "increasing"
                directionLabels(i) = "↓↑";
            elseif chronic == "no change"
                directionLabels(i) = "↓x";
            elseif chronic == "decreasing"
                directionLabels(i) = "↓↓";
            end
        end
    end

    directions = ["↑↑", "↑x", "↑↓", "x↑", "xx", "x↓", "↓↑", "↓x", "↓↓"];
    numDirections = numel(directions);

    % Unique cell types
    [uniqueTypes, ~, typeIdx] = unique(cellTypes);
    numTypes = numel(uniqueTypes);
    
    % Initialize count matrix: rows = cell types, columns = direction labels
    directionCounts = zeros(numTypes, numDirections);
    clusterID = size(cellTypes,1);

    % Count occurrences
    for i = 1:length(directionLabels)
        row = typeIdx(i);  % cell type index
        col = find(directions == directionLabels(i));  % direction label index
        directionCounts(row, col) = directionCounts(row, col) + 1;
        clusterID(i,1) = col;
    end

end

%     Group	pre→during	during→post
%     1	 ↑	↑
%     2	 ↑	x
%     3	 ↑	↓
%     4	 x	↑
%     5	 x	x
%     6	 x	↓
%     7	 ↓	↑
%     8	 ↓	x
%     9	 ↓	↓


