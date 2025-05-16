function sortedFileNames = sortFileNames(fileNames)
% sortFileNamesByAnimalDateFOV sorts filenames by animal ID, date, and FOV number.
%
% INPUT:
%   fileNames - cell array of strings (e.g., {'Pf4Ai162-2_221130_FOV6', ...})
%
% OUTPUT:
%   sortedFileNames - cell array of sorted filenames

    uniqueFiles = unique(fileNames);  % Ensure uniqueness
    sortKey = zeros(length(uniqueFiles), 3);  % Preallocate sort key matrix

    for file = 1:length(uniqueFiles)
        filename = uniqueFiles{file};

        % Extract number after 'Pf4Ai162-'
        numberAfterPrefix = sscanf(filename, 'Pf4Ai162-%d', 1);

        % Extract 6-digit date
        dateStr = regexp(filename, '\d{6}', 'match', 'once');
        dateNum = str2double(dateStr);

        % Extract FOV number
        fovMatch = regexp(filename, 'FOV(\d+)', 'tokens', 'once');
        if ~isempty(fovMatch)
            fovNumber = str2double(fovMatch{1});
        else
            fovNumber = NaN;  % Handle missing FOV gracefully
        end

        sortKey(file, :) = [numberAfterPrefix, dateNum, fovNumber];
    end

    % Sort based on the composite key
    [~, idx] = sortrows(sortKey);
    sortedFileNames = uniqueFiles(idx);
end
