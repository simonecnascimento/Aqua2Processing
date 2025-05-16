function plotClusterHeatmap_FOVnormalized(combinedTable, includedClusters, sortedFileNames, savePath)
% plotClusterHeatmap - Plots a heatmap of calcium traces grouped by cluster,
%                      with normalization by FOV and optional z-scoring.
%
% Inputs:
%   combinedTable    - Table with cell info including dFF traces, cluster IDs, and FOV names
%   includedClusters - Vector of cluster IDs to include (e.g. [2,4,8])
%   sortedFileNames  - Cell array of FOV names in desired order
%   savePath         - Optional folder to save figure

    % === Parameters ===
    baselineRange = 1:1854;
    traceCol = combinedTable{:,"dFF"};
    fileNames = combinedTable{:,"fileNameColumn"};
    clusterIDs = combinedTable{:,"clusterID"};
    normalizedTraces = cell(height(combinedTable), 1);

    % === Normalize each FOV by its baseline ===
    for i = 1:numel(sortedFileNames)
        thisFOV = sortedFileNames{i};
        fovIdx = strcmp(fileNames, thisFOV);
        traces = traceCol(fovIdx);

        % Get FOV-wide baseline from mean across all cells
        allBaselines = zeros(numel(traces), 1);
        for j = 1:numel(traces)
            allBaselines(j) = mean(traces{j}(baselineRange));
        end
        fovBaseline = mean(allBaselines);

        % Normalize traces by FOV baseline
        fovIndices = find(fovIdx);
        for j = 1:numel(traces)
            normalizedTraces{fovIndices(j)} = (traces{j} - fovBaseline) / fovBaseline;
        end
    end

    % === Collect traces by cluster ===
    selectedData = {};
    sortedClusters = [];
    for i = 1:length(includedClusters)
        clust = includedClusters(i);
        idx = find(clusterIDs == clust);
        selectedData = [selectedData; normalizedTraces(idx)];
        sortedClusters = [sortedClusters; repmat(clust, length(idx), 1)];
    end

    % === Convert to matrix and apply z-score (per row/cell) ===
    dFF_Matrix = cell2mat(selectedData);  % N x T
  
    figure;
    imagesc(dFF_Matrix);
    colormap(flipud(gray));
    colorbar;
    xlabel('Time (frames)');
    ylabel('Cells (sorted by cluster)');
    title(['Z-scored Heatmap (clusters: ', num2str(includedClusters), ')']);
    hold on;

    % Add cluster boundaries
    clusterBoundaries = find(diff(sortedClusters));
    for i = 1:length(clusterBoundaries)
        yline(clusterBoundaries(i)+0.5, 'Color', 'k', 'LineWidth', 2);
    end

    % Time marker (optional)
    xline(1854, '--r', 'LineWidth', 1.5);

    % Save if needed
    if nargin >= 4 && ~isempty(savePath)
        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end
        filename = fullfile(savePath, 'clusterHeatmap_zscore.eps');
        saveas(gcf, filename);
        fprintf('Figure saved to: %s\n', filename);
    end
end

