function plotClusterDistributionByFOV(combinedTable_clusters, savePath)
    % This function plots pie charts for each FOV, showing the distribution 
    % of clusters. It saves the figure to the specified output file.
    %
    % Parameters:
    %   combinedTable_clusters: A table with columns 'fileNameColumn' and 'clusterID'.
    %   outputFileName: The name (and path) of the output file where the figure should be saved.
    
%     % Get unique filenames
%     uniqueFiles = unique(combinedTable_clusters.fileNameColumn);
% 
%     % Initialize an array to store the extracted parts for sorting
%     sortKey = [];
% 
%     % Loop through all filenames to extract the parts for sorting
%     for file = 1:length(uniqueFiles)
%         filename = uniqueFiles{file};
%         
%         % Extract the number after "Pf4Ai162-" (e.g., '2' from 'Pf4Ai162-2')
%         numberAfterPrefix = sscanf(filename, 'Pf4Ai162-%d', 1);
%         
%         % Extract the date (e.g., '221130' from 'Pf4Ai162-2_221130_FOV6')
%         dateStr = regexp(filename, '\d{6}', 'match', 'once');
%         
%         % Extract the FOV number (e.g., 'FOV6' from 'Pf4Ai162-2_221130_FOV6')
%         fovNumber = sscanf(filename, 'Pf4Ai162-%*d_%*d_FOV%d', 1);
%         
%         % Store the extracted values in a matrix for sorting
%         sortKey = [sortKey; numberAfterPrefix, str2double(dateStr), fovNumber];
%     end
% 
%     % Sort by the three columns: numberAfterPrefix, date, fovNumber
%     [~, idx] = sortrows(sortKey);
%     sortedFileNames = uniqueFiles(idx);

    fileNames = combinedTable_clusters.fileNameColumn;
    sortedFileNames = sortFileNames(fileNames);

    % Create FOV labels (e.g., 'FOV 1', 'FOV 2', ..., 'FOV 13')
    FOVlabels = arrayfun(@(x) ['FOV ' num2str(x)], 1:length(sortedFileNames), 'UniformOutput', false);

    % Create containers.Map as dictionary
    FOVmap = containers.Map(sortedFileNames, FOVlabels);

    % Get unique FOVs (just numbers)
    FOVnumbers = 1:length(sortedFileNames);

    % Set up figure layout for up to 13 FOVs
    fig = figure;
    tiledlayout(4, 4); % 4x4 grid of plots (adjustable if more than 13 FOVs)

    % Loop through each FOV
    for i = 1:length(FOVnumbers)
        FOVID = FOVnumbers(i);
        currentFile = sortedFileNames{FOVID};
        % Filter rows that belong to this FOV
        fovData = combinedTable_clusters(contains(combinedTable_clusters.fileNameColumn, currentFile), :);
        
        % Get the cluster IDs for this FOV
        clusterIDs = fovData.eventRate_clusterID;

        % Count occurrences of each cluster for this FOV
        [uniqueClusters, ~, idx] = unique(clusterIDs);
        counts = accumarray(idx, 1);

        % Plot pie chart
        % Remove clusters with count == 0
        validIdx = counts > 0;
        counts = counts(validIdx);
        clusterLabels = arrayfun(@(x) [num2str(x)], uniqueClusters(validIdx), 'UniformOutput', false);
        
        % Plot pie chart
        nexttile;
        h = pie(counts, clusterLabels);

        % Set font size for labels and title
        for j = 1:length(h)
            if isprop(h(j), 'Text')
                h(j).FontSize = 8; % Adjust the font size as needed
            end
        end
        title(['FOV ' num2str(FOVID)], 'FontSize', 10); % Adjust title font size
    end

    if nargin >= 2 && ~isempty(savePath)
        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end
        filename = fullfile(savePath, 'Cluster_Distribution_FOV.png');
        saveas(fig, filename);
        fprintf('Figure saved to: %s\n', filename);
    end
end
