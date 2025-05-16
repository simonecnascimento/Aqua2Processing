function sortedFileNames = sortFileNamesPf4Ai162(fileNames)
    % sortFileNamesPf4Ai162 sorts filenames based on numeric components
    %
    % Input:
    %   fileNames - cell array of filenames (e.g., {'Pf4Ai162-2_221130_FOV6', ...})
    %
    % Output:
    %   sortedFileNames - filenames sorted by animal number, date, and FOV

    % Initialize array to store extracted parts for sorting
    sortKey = [];

    for file = 1:length(fileNames)
        filename = fileNames{file};

        % Extract the number after "Pf4Ai162-"
        numberAfterPrefix = sscanf(filename, 'Pf4Ai162-%d', 1);

        % Extract the 6-digit date
        dateStr = regexp(filename, '\d{6}', 'match', 'once');
        dateNum = str2double(dateStr);

        % Extract FOV number
        fovNumber = sscanf(filename, 'Pf4Ai162-%*d_%*d_FOV%d', 1);

        % Append to sortKey
        sortKey = [sortKey; numberAfterPrefix, dateNum, fovNumber];
    end

    % Sort by the composite key
    [~, idx] = sortrows(sortKey);
    sortedFileNames = fileNames(idx);
end
