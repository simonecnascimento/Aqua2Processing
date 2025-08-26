function startingFrames_byCell = extractStartingFramesByCell(eventsByCell_experiment, startingFrame)
% Extracts starting frames per cell and per condition (preCSD, duringCSD, postCSD, baseline)
%
% Inputs:
%   - eventsByCell_experiment: cell array with events per cell, columns 2–5 are event lists
%   - startingFrame: a numeric vector or array mapping event indices to frame numbers
%
% Output:
%   - startingFrames_byCell: cell array (numCells x 4) with starting frames per condition

startingFrames_byCell = cell(size(eventsByCell_experiment, 1), 5);  % Add one extra column for cell ID

for i = 1:size(eventsByCell_experiment, 1)
    % Copy cell ID or label into first column
    startingFrames_byCell{i, 1} = eventsByCell_experiment{i, 1};

    for j = 1:4  % Process columns 2 to 5 of eventsByCell_experiment
        eventList = eventsByCell_experiment{i, j+1};

        if isempty(eventList)
            startingFrames_byCell{i, j+1} = [];

        elseif isnumeric(eventList)
            % Assume eventList contains event indices; map to starting frames
            try
                startingFrames_byCell{i, j+1} = startingFrame(eventList);
            catch
                warning('Index exceeds startingFrame bounds at cell (%d, %d)', i, j+1);
                startingFrames_byCell{i, j+1} = [];
            end

        else
            warning('Unexpected format in eventsByCell_experiment at cell (%d, %d)', i, j+1);
            startingFrames_byCell{i, j+1} = [];
        end
    end
end

