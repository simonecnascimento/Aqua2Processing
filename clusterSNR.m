function [idx, clustered_cells] = clusterSNR(snr_values, maxClusters, numClusters)
    % CLUSTERSNR Clusters cells based on SNR values using k-means
    %
    %   [idx, clustered_cells] = clusterSNR(snr_values, maxClusters, numClusters)
    %
    %   Inputs:
    %       snr_values  - Column vector of SNR values (one per cell)
    %       maxClusters - Maximum number of clusters to test for Elbow method
    %       numClusters - Number of clusters to use for final k-means clustering
    %
    %   Outputs:
    %       idx            - Cluster indices for each cell
    %       clustered_cells - Cell array containing indices of cells in each cluster

    % Step 1: Determine optimal number of clusters using Elbow method
    sum_of_squares = zeros(maxClusters, 1);
    for k = 1:maxClusters
        [~, ~, sumd] = kmeans(snr_values, k, 'Replicates', 5);
        sum_of_squares(k) = sum(sumd);
    end
    
    % Plot the Elbow curve
    figure;
    plot(1:maxClusters, sum_of_squares, '-o', 'LineWidth', 2);
    xlabel('Number of Clusters'); ylabel('Sum of Squared Distances');
    title('Elbow Method for Optimal Clusters');
    grid on;
    
    % Step 2: Apply k-means clustering with specified numClusters
    [idx, ~] = kmeans(snr_values, numClusters, 'Replicates', 5);
    
    % Assign cells to clusters
    clustered_cells = cell(numClusters, 1);
    for k = 1:numClusters
        clustered_cells{k} = find(idx == k);
        fprintf('Cluster %d: %d cells\n', k, length(clustered_cells{k}));
    end
    
    % Plot histogram of SNR values
    figure;
    histogram(snr_values, 30);
    xlabel('SNR'); ylabel('Number of Cells');
    title(sprintf('SNR Distribution (Clusters = %d)', numClusters));
end
