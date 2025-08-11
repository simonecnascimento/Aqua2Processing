function plotFOVDistributionByCluster(combinedTable_clusters, savePath)
    % This function plots pie charts for each cluster, showing the distribution 
    % of different FOVs. It saves the figure to the specified output file.
    %
    % Parameters:
    %   combinedTable_clusters: A table with columns 'fileNameColumn' and 'clusterID'.
    %   outputFileName: The name (and path) of the output file where the figure should be saved.
    
    fileNames = combinedTable_clusters.fileNameColumn;
    sortedFileNames = sortFileNames(fileNames);

    % Create FOV labels (e.g., 'FOV 1', 'FOV 2', ..., 'FOV 13')
    FOVlabels = arrayfun(@(x) ['FOV ' num2str(x)], 1:length(sortedFileNames), 'UniformOutput', false);

    % Create containers.Map as dictionary
    FOVmap = containers.Map(sortedFileNames, FOVlabels);

    % Get unique clusters
    uniqueClusters = unique(combinedTable_clusters.eventRate_clusterID);

    % Set up figure layout for up to 9 clusters
    fig = figure;
    tiledlayout(3,3); % 3x3 grid of plots

    % Loop through each cluster
    for i = 1:length(uniqueClusters)
        clusterID = uniqueClusters(i);

        % Filter rows that belong to this cluster
        clusterData = combinedTable_clusters(combinedTable_clusters.eventRate_clusterID == clusterID, :);

        rawNames = clusterData.fileNameColumn;

        % Map each filename to its FOV label using the dictionary
        labelsMapped = cellfun(@(x) FOVmap(x), rawNames, 'UniformOutput', false);

        % Count occurrences per FOV
        [uniqueLabels, ~, idx] = unique(labelsMapped);
        counts = accumarray(idx, 1);

        validIdx = counts > 0;
        counts = counts(validIdx);
        FOV_Labels = uniqueLabels(validIdx);  % no need for num2str

        % Plot pie chart
        nexttile;
        h = pie(counts, uniqueLabels);

        % Set font size for labels and title
        for j = 1:length(h)
            % Adjust font size of pie chart labels
            if isprop(h(j), 'Text')
                h(j).FontSize = 8; % adjust the font size as needed
            end
        end
        title(['Cluster ' num2str(clusterID)], 'FontSize', 10); % adjust title font size
    end

    % Save figure if savePath is provided
    if nargin >= 2 && ~isempty(savePath)
        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end
        filename = fullfile(savePath, 'FOV_Distribution_Cluster.png');
        saveas(fig, filename);
        fprintf('Figure saved to: %s\n', filename);
    end

%     % Save figure as the specified output file
%     saveas(gcf, outputFileName);  % You can change the file path and extension (e.g., .jpg, .fig)
end
