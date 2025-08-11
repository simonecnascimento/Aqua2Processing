function [events_byCell_all, data_byCell_all, data_byFOV_all] = calculateEventRate(paramTables_all_byCell, paramTables_all_byFOV, durationCSDmin)
    % CALCULATEEVENTRATE Calculate event rates for cells and FOVs.
    %
    % Inputs:
    %   paramTables_all_byCell - Table with event counts by cell
    %   paramTables_all_byFOV - Table with event counts by FOV
    %   durationMinutes - Vector specifying the duration (in minutes) for each column
    %
    % Outputs:
    %   eventHz_byCell - Event rates by cell (Hz)
    %   eventHz_byFOV - Event rates by FOV (Hz)

    % Convert duration to seconds
    durationMinutes = [30, durationCSDmin, 30]; % Duration in minutes for each column, column1=cellID
    durationSeconds = durationMinutes * 60;

    % empty structures
    data_byCell_all = [];
    data_byFOV_all = [];

    %
    events_byCell_all = [];

    % byCell
    for e = 1:length(paramTables_all_byCell)

        % byCell Hz
        data_byCell = table2cell(paramTables_all_byCell(e).NumberOfEvents(:,1:4));
        emptyEvents_byCell = cellfun(@isempty, data_byCell);       % Find empty values
        data_byCell(emptyEvents_byCell) = {0};                   % Replace with 0
        data_byCell = cell2mat(data_byCell);          % Convert to numeric matrix
        
        % raw event count
        rowToRemove = (data_byCell(:, 2) == 0 & data_byCell(:, 3) == 0 & data_byCell(:, 4) == 0);
        data_byCell(rowToRemove, :) = [];
        events_byCell_all = [events_byCell_all; data_byCell];
        
        % byCell event
        eventHz_byCell_phases = data_byCell(:,2:4) ./ durationSeconds;         % Calculate event rate (Hz)
        eventHz_byCell = [data_byCell(:,1), eventHz_byCell_phases]; % add cellID column

        rowsToRemove = (eventHz_byCell(:, 2) == 0 & eventHz_byCell(:, 3) == 0 & eventHz_byCell(:, 4) == 0);
        eventHz_byCell(rowsToRemove, :) = [];

        data_byCell_all = [data_byCell_all; eventHz_byCell];

    end

    % Logical conditions for classification of cells
    risingEvents = data_byCell_all(:, 3) > 2 * data_byCell_all(:, 2);
    decreasingEvents = data_byCell_all(:, 3) < 0.5 * data_byCell_all(:, 2);
    noChangeEvents = ~(risingEvents | decreasingEvents);
    
    % Initialize a column for trends as a cell array to store strings
    eventTrends = strings(size(data_byCell_all, 1), 1);
    
    % Assign trends using logical indexing
    eventTrends(risingEvents) = "rising";
    eventTrends(decreasingEvents) = "decreasing";
    eventTrends(noChangeEvents) = "no change";
    
    % Add the trend column to the matrix 
    data_byCell_all = [data_byCell_all, eventTrends];
    
    increasingCells = data_byCell_all(risingEvents, :);
    decreasingCells = data_byCell_all(decreasingEvents, :);
    noChangeCells = data_byCell_all(noChangeEvents, :);

    % check how many cells have no recovered
    noRecoveryCells = decreasingCells(decreasingCells(:,4) == "0", :);
    recoveryCells = decreasingCells(decreasingCells(:,4) ~= "0", :);
    
    % byFOV
    data_byFOV = table2cell(paramTables_all_byFOV(:,1:3));
    emptyEvents_byFOV = cellfun(@isnan, data_byFOV);        % Find NaN values
    data_byFOV(emptyEvents_byFOV) = {0};                    % Replace NaN with 0
    data_byFOV = cell2mat(data_byFOV);           % Convert to numeric matrix
    eventHz_byFOV = data_byFOV ./ durationSeconds;          % Calculate event rate (Hz)
    data_byFOV_all = eventHz_byFOV;

    % Logical conditions for classification of events
    risingFOV = data_byFOV_all(:, 2) > 2 * data_byFOV_all(:, 1);
    decreasingFOV = data_byFOV_all(:, 2) < 0.5 * data_byFOV_all(:, 1);
    noChangeFOV = ~(risingFOV | decreasingFOV);
    
    % Initialize a column for trends as a cell array to store strings
    eventTrendsFOV = strings(size(data_byFOV_all, 1), 1);
    
    % Assign trends using logical indexing
    eventTrendsFOV(risingFOV) = "rising";
    eventTrendsFOV(decreasingFOV) = "decreasing";
    eventTrendsFOV(noChangeFOV) = "no change";
    
    % Add the trend column to the matrix 
    data_byFOV_all = [data_byFOV_all, eventTrendsFOV];
    
    increasingFOV = data_byFOV_all(risingFOV, :);
    decreasingFOV = data_byFOV_all(decreasingFOV, :);
    noChangeFOV = data_byFOV_all(noChangeFOV, :);
end
