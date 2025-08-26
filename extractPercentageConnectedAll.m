function percentageConnectedAll_flat = extractPercentageConnectedAll(NetworkTable)
% extractPercentageConnectedAll - Extracts and concatenates all numeric values
% from NetworkTable.percentageConnectedAll for distribution analysis.
%
% INPUT:
%   NetworkTable - struct containing field:
%                  .percentageConnectedAll : cell array of Nx1 doubles or empty arrays
%
% OUTPUT:
%   percentageConnectedAll_flat - combined numeric column vector of all values
%
% EXAMPLE USAGE:
%   percentageConnectedAll_flat = extractPercentageConnectedAll(NetworkTable);
%   histogram(percentageConnectedAll_flat);

    % Ensure the field exists
    if ~isfield(NetworkTable, 'percentageConnectedAll')
        error('NetworkTable does not contain the field "percentageConnectedAll".');
    end

    % Extract the cell array
    pcAll = NetworkTable.percentageConnectedAll;

    % Filter out empty cells and concatenate non-empty numeric arrays
    percentageConnectedAll_flat = vertcat(pcAll{~cellfun(@isempty, pcAll)});

end
