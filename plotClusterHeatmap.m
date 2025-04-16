function plotClusterHeatmap(combinedTable, clusterCol, traceCol, includedClusters)
% plotClusterHeatmap - Plots a heatmap of calcium traces grouped by cluster
%
% Syntax:
%   plotClusterHeatmap(combinedTable, clusterCol, traceCol, includedClusters)
%
% Inputs:
%   combinedTable      - MATLAB table containing all cell data
%   clusterCol         - Column index of the cluster IDs
%   traceCol           - Column index of the calcium traces (1xT double arrays)
%   includedClusters   - Vector of cluster IDs to include (e.g. [2,4,8,9])
%
% Example:
%   plotClusterHeatmap(myTable, 17, 15, [2, 4, 8, 9]);

    % Initialize
    selectedData = {};
    selectedClusters = [];

    % Loop through each cluster
    for i = 1:length(includedClusters)
        clust = includedClusters(i);
        idx = find(combinedTable{:, clusterCol} == clust);

        % Limit to max 10 cells per cluster
        if numel(idx) >  20
            idx = idx(1:20); %  idx = idx(1:10); or use randperm(numel(idx), 10) for random
        end

        % Append
        selectedData = [selectedData; combinedTable{idx, traceCol}];
        selectedClusters = [selectedClusters; repmat(clust, length(idx), 1)];
    end

    % Convert to matrix
    dFF_Matrix = cell2mat(selectedData);

    % Sort by cluster
    [sortedClusters, sortIdx] = sort(selectedClusters);
    sorted_dFF = dFF_Matrix(sortIdx, :);

    % Plot
    figure;
    imagesc(sorted_dFF);
    colormap(flipud(gray));
    colorbar;
    xlabel('Time (frames)');
    ylabel('Cells (sorted by cluster)');
    title(['Heatmap (max 10 cells per cluster): ', num2str(includedClusters)]);

    % Cluster dividers
    hold on;
    clusterBoundaries = find(diff(sortedClusters));
    for i = 1:length(clusterBoundaries)
        yline(clusterBoundaries(i)+0.5, 'Color', 'k', 'LineWidth', 1);
    end
end

%%

%     % Extract cluster IDs
%     clusterIDs = combinedTable{:, clusterCol};
% 
%     % Find selected clusters
%     selectedIdx = ismember(clusterIDs, includedClusters);
% 
%     % Extract traces and cluster IDs
%     filteredData = combinedTable{selectedIdx, traceCol};
%     filteredClusters = clusterIDs(selectedIdx);
% 
%     % Convert to matrix
%     dFF_Matrix = cell2mat(filteredData);
% 
%     % Sort rows by cluster ID
%     [sortedClusters, sortIdx] = sort(filteredClusters);
%     sorted_dFF = dFF_Matrix(sortIdx, :);
% 
%     % Plot heatmap
%     figure;
%     imagesc(sorted_dFF);
%     colormap(flipud(gray));
%     colorbar;
%     xlabel('Time (frames)');
%     ylabel('Cells (sorted by cluster)');
%     title(['Heatmap of Clusters: ', num2str(includedClusters)]);
% 
%     % Add cluster boundaries
%     hold on;
%     clusterBoundaries = find(diff(sortedClusters));
%     for i = 1:length(clusterBoundaries)
%         yline(clusterBoundaries(i) + 0.5, 'Color', 'b', 'LineWidth', 2);
%     end

