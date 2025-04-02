function cellTable = categorize_waves(combinedTable_NM)

    % Calculate quartiles for Duration
    medianDuration = median(combinedTable_NM.("Duration 10% to 10%"));
    Q1_duration = prctile(combinedTable_NM.("Duration 10% to 10%"), 25);
    Q3_duration = prctile(combinedTable_NM.("Duration 10% to 10%"), 75);

    % Identify indices for duration categories
    highDurationIdx = find(combinedTable_NM.("Duration 10% to 10%") >= medianDuration & ...
                           combinedTable_NM.("Duration 10% to 10%") <= Q3_duration);
    lowDurationIdx = find(combinedTable_NM.("Duration 10% to 10%") >= Q1_duration & ...
                          combinedTable_NM.("Duration 10% to 10%") < medianDuration);

    % Calculate quartiles for Frequency
    medianFrequency = median(combinedTable_NM.("Number of Events"));
    Q1_frequency = prctile(combinedTable_NM.("Number of Events"), 25);
    Q3_frequency = prctile(combinedTable_NM.("Number of Events"), 75);

    % Identify indices for frequency categories
    highFrequencyIdx = find(combinedTable_NM.("Number of Events") >= medianFrequency & ...
                            combinedTable_NM.("Number of Events") <= Q3_frequency);
    lowFrequencyIdx = find(combinedTable_NM.("Number of Events") >= Q1_frequency & ...
                           combinedTable_NM.("Number of Events") < medianFrequency);

    % Calculate quartiles for max dFF
    medianMaxDFF = median(combinedTable_NM.("Max dFF"));
    Q1_maxdFF = prctile(combinedTable_NM.("Max dFF"), 25);
    Q3_maxdFF = prctile(combinedTable_NM.("Max dFF"), 75);

    % Identify indices for max dFF categories
    highMaxDFFIdx = find(combinedTable_NM.("Max dFF") >= medianMaxDFF & ...
                         combinedTable_NM.("Max dFF") <= Q3_maxdFF);
    lowMaxDFFIdx = find(combinedTable_NM.("Max dFF") >= Q1_maxdFF & ...
                        combinedTable_NM.("Max dFF") < medianMaxDFF);

    % Create an empty table with 503 rows
    numCells = size(combinedTable_NM, 1);
    cellTable = table((1:numCells)', nan(numCells,1), nan(numCells,1), nan(numCells,1), nan(numCells,1), ...
        'VariableNames', {'CellNumber', 'CellType', 'Duration', 'Frequency', 'Max_dFF'});

    % Assign categories to the table
    for x = 1:numCells

        cellTable.CellType(x) = combinedTable_NM.("Cell location (0,perivascular;1,adjacent;2,none)")(x,1);

        % Duration Category
        if ismember(x, highDurationIdx)
            cellTable.Duration(x) = 2;  % High Duration
        elseif ismember(x, lowDurationIdx)
            cellTable.Duration(x) = 1;  % Low Duration
        else
            cellTable.Duration(x) = 0;  % Other
        end

        % Frequency Category
        if ismember(x, highFrequencyIdx)
            cellTable.Frequency(x) = 2;  % High Frequency
        elseif ismember(x, lowFrequencyIdx)
            cellTable.Frequency(x) = 1;  % Low Frequency
        else
            cellTable.Frequency(x) = 0;  % Other
        end

        % Max dFF Category
        if ismember(x, highMaxDFFIdx)
            cellTable.Max_dFF(x) = 2;  % High max dFF
        elseif ismember(x, lowMaxDFFIdx)
            cellTable.Max_dFF(x) = 1;  % Low max dFF
        else
            cellTable.Max_dFF(x) = 0;  % Other
        end
    end
end
