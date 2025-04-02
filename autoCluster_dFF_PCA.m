function [idx, clustered_cells, numClusters] = autoCluster_dFF_PCA(dFF_all, numPCs, maxClusters)
    % AUTOCLUSTER_DFF_PCA Clusters dFF signals using PCA and auto-selects clusters
    %
    % Inputs:
    %   dFF_all     - Matrix of dFF signals (cells x time points)
    %   numPCs      - Number of principal components to use
    %   maxClusters - Maximum number of clusters to test (for Elbow method)
    %
    % Outputs:
    %   idx             - Cluster indices for each cell
    %   clustered_cells - Cell array of indices for each cluster
    %   numClusters     - Optimal number of clusters
    
    % Perform PCA
    [coeff, score, ~] = pca(dFF_all);  
    
    % Reduce data to selected PCs
    reduced_data = score(:, 1:numPCs);
    
    % Find optimal number of clusters using the Elbow method
    sum_of_squares = zeros(maxClusters, 1);
    for k = 1:maxClusters
        [~, ~, sumd] = kmeans(reduced_data, k, 'Replicates', 5);
        sum_of_squares(k) = sum(sumd);
    end

    % Determine the "elbow" point
    diff_sos = diff(sum_of_squares);      % First derivative
    diff_sos2 = diff(diff_sos);           % Second derivative
    [~, numClusters] = max(abs(diff_sos2)); % Elbow at max curvature

    % Apply k-Means clustering with optimal numClusters
    [idx, ~] = kmeans(reduced_data, numClusters, 'Replicates', 5);
    
    % Assign clusters
    clustered_cells = cell(numClusters, 1);
    for k = 1:numClusters
        clustered_cells{k} = find(idx == k);
    end
    
    % Plot Elbow method
    figure;
    plot(1:maxClusters, sum_of_squares, '-o', 'LineWidth', 2);
    hold on;
    plot(numClusters, sum_of_squares(numClusters), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    xlabel('Number of Clusters'); ylabel('Sum of Squared Distances');
    title('Elbow Method for Optimal Clusters');
    grid on;

    % Plot mean waveforms for each cluster
    figure;
    hold on;
    colors = lines(numClusters);
    for k = 1:numClusters
        plot(mean(dFF_all(clustered_cells{k}, :), 1), 'Color', colors(k, :), 'LineWidth', 2);
    end
    title('Clustered dFF Signals');
    xlabel('Time'); ylabel('dF/F');
    legend(arrayfun(@(x) sprintf('Cluster %d', x), 1:numClusters, 'UniformOutput', false));
    hold off;
end
