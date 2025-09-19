function event_cell_Table = fillMissingCellRows(event_cell_Table, data_CFU)
% fillMissingCellRows ensures that event_cell_Table has rows for all cell IDs
%
% Input:
%   event_cell_Table - cell array {cellID, eventList, signal}
%
% Output:
%   event_cell_Table - updated table with empty rows for missing cellIDs

    % Find maximum cell number
    %maxCell = max(cell2mat(event_cell_Table(:,1)));
    maxCell = numel(data_CFU.cfuInfo1(:,1));

    % Preallocate full table
    newTable = cell(maxCell, size(event_cell_Table,2));

    % Fill with defaults
    for i = 1:maxCell
        newTable{i,1} = i;   % cell ID
        newTable{i,2} = [];  % empty event list
        newTable{i,3} = [];  % empty signal
    end

    % Insert existing data
    for r = 1:size(event_cell_Table,1)
        id = event_cell_Table{r,1};
        newTable{id,2} = event_cell_Table{r,2};
        newTable{id,3} = event_cell_Table{r,3};
    end

    % Return updated table
    event_cell_Table = newTable;
end
