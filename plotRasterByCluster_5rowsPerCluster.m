function plotRasterByCluster_5rowsPerCluster(rasterTable)
    % Plots raster with all events (preCSD, duringCSD, postCSD) in black
    % and arranges cells so that each cluster occupies 5 rows (padding empty rows if needed).
    % X axis limited to 3772.
    
    clusters = unique(rasterTable.clusterID);
    numClusters = length(clusters);
    
    rowsPerCluster = 5;
    totalRows = rowsPerCluster * numClusters;
    
    figure;
    hold on;
    
    rowPositions = nan(height(rasterTable), 1); % To store y positions for each cell
    
    for cIdx = 1:numClusters
        clusterID = clusters(cIdx);
    
        % Find rows for this cluster
        clusterRows = find(rasterTable.clusterID == clusterID);
    
        % Assign rows within cluster block:
        % If fewer than 5 cells, remaining rows stay empty
        for i = 1:length(clusterRows)
            % y position = cluster start row + offset (0 to 4)
            rowPositions(clusterRows(i)) = (cIdx - 1) * rowsPerCluster + i;
        end
    end
    
    % Now plot events per cell at assigned row positions
    for i = 1:height(rasterTable)
        times = [rasterTable.preCSD{i}, rasterTable.duringCSD{i}, rasterTable.postCSD{i}];
        y = rowPositions(i);
        if isnan(y)
            % Cell not assigned a row (should not happen)
            continue;
        end
        for t = times
            if t <= 3772
                line([t t], [y-0.4, y+0.4], 'Color', 'k', 'LineWidth', 1);
            end
        end
    end
    
       % Vertical lines at frames 1854 and 1855
        xline(1854, '--r', 'LineWidth', 1.5);

    xlabel('Frame');
    ylabel('Cluster (5 rows each)');
    xlim([0 3772]);
    ylim([0 totalRows+1]);
    title('Raster Plot: 5 rows per Cluster');
    set(gca, 'YDir', 'reverse');
    box on;
    hold off;
end
