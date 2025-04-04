function [idx, clustered_cells, optimalK, count_k1, count_k2] = clusterSNR(snr_values, maxClusters, dFF_all)
    % CLUSTERSNR Clusters cells based on precomputed SNR values using k-means
    %
    %   [idx, clustered_cells, optimalK] = clusterSNR(snr_values, maxClusters)
    %
    %   Inputs:
    %       snr_values  - Column vector of SNR values (one per cell)
    %       maxClusters - Maximum number of clusters to test for Elbow method
    %
    %   Outputs:
    %       idx            - Cluster indices for each cell
    %       clustered_cells - Cell array containing indices of cells in each cluster
    %       optimalK       - Optimal number of clusters based on the Elbow Method
    
    % Step 1: Compute within-cluster sum of squares for different k values (Elbow Method)
    sum_of_squares = zeros(maxClusters, 1);
    for k = 1:maxClusters
        [~, ~, sumd] = kmeans(snr_values, k, 'Replicates', 5);
        sum_of_squares(k) = sum(sumd);
    end
    
    % Step 2: Find optimal number of clusters using second derivative (Elbow Method)
    first_diffs = diff(sum_of_squares);
    second_diffs = diff(first_diffs);

    % Check if second differences have enough data and make sure the difference is significant
    if length(second_diffs) > 1
        [~, optimalK] = max(abs(second_diffs)); % +1 because diff reduces size
        optimalK = optimalK + 1; 
    else
        % If there’s not enough second difference data, fall back on a simple heuristic
        [~, optimalK] = min(first_diffs);
        optimalK = optimalK + 1; % Adjust to the correct number of clusters
    end
      
    % Step 3: Plot the Elbow curve
    figure;
    plot(1:maxClusters, sum_of_squares, '-o', 'LineWidth', 2);
    hold on;
    plot(optimalK, sum_of_squares(optimalK), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Number of Clusters'); ylabel('Sum of Squared Distances');
    title(sprintf('Elbow Method: Optimal Clusters = %d', optimalK));
    grid on;
    hold off;
    
    % Step 4: Apply k-means clustering with the chosen number of clusters
    [idx, ~] = kmeans(snr_values, optimalK, 'Replicates', 5, 'MaxIter', 300, 'Start', 'plus');
    count_k1 = sum(idx == 1);
    count_k2 = sum(idx == 2);
    
    % Step 5: Assign cells to clusters
    clustered_cells = cell(optimalK, 1);
    for k = 1:optimalK
        clustered_cells{k} = find(idx == k);
        fprintf('Cluster %d: %d cells\n', k, length(clustered_cells{k}));
    end
    
    % Step 6: Plot histogram of SNR values with cluster assignments
    figure;
    histogram(snr_values, 70);
    xlabel('SNR'); ylabel('Number of Cells');
    title(sprintf('SNR Distribution (Optimal Clusters = %d)', optimalK));

    % Find the global min and max
    minData = min(dFF_all(:));  % Global min
    maxData = max(dFF_all(:));  % Global max

    % Plot waveforms by cluster
    figure;
    for k = 1:optimalK
        subplot(optimalK, 1, k);
        plot(dFF_all(idx == k, :)', 'LineWidth', 1);
        title(sprintf('Cluster %d', k));
        xlabel('Time'); ylabel('\DeltaF/F');
        ylim([minData, maxData]);  
        xlim([0,900]);
        box off;
    end

    % Plot heatmap by cluster
    for k = 1:optimalK
        data = dFF_all(idx == k, :);
        data = data(1:count_k2, :); 
         % Use imagesc instead of heatmap
        figure;
        imagesc(data); % or imagesc(data) if normalization is not needed
        caxis([minData, maxData]); % Normalize the color scale
        colormap parula;
        colorbar;
        title(sprintf('Cluster %d Heatmap dFF', k));
        xlabel('Time');
        ylabel('Cell Index');
    end

end



% function [idx, clustered_cells] = clusterSNR(snr_values, maxClusters, numClusters)
%     % CLUSTERSNR Clusters cells based on SNR values using k-means
%     %
%     %   [idx, clustered_cells] = clusterSNR(snr_values, maxClusters, numClusters)
%     %
%     %   Inputs:
%     %       snr_values  - Column vector of SNR values (one per cell)
%     %       maxClusters - Maximum number of clusters to test for Elbow method
%     %       numClusters - Number of clusters to use for final k-means clustering
%     %
%     %   Outputs:
%     %       idx            - Cluster indices for each cell
%     %       clustered_cells - Cell array containing indices of cells in each cluster
% 
%     % Step 1: Determine optimal number of clusters using Elbow method
%     sum_of_squares = zeros(maxClusters, 1);
%     for k = 1:maxClusters
%         [~, ~, sumd] = kmeans(snr_values, k, 'Replicates', 5);
%         sum_of_squares(k) = sum(sumd);
%     end
%     
%     % Plot the Elbow curve
%     figure;
%     plot(1:maxClusters, sum_of_squares, '-o', 'LineWidth', 2);
%     xlabel('Number of Clusters'); ylabel('Sum of Squared Distances');
%     title('Elbow Method for Optimal Clusters');
%     grid on;
%     
%     % Step 2: Apply k-means clustering with specified numClusters
%     [idx, ~] = kmeans(snr_values, numClusters, 'Replicates', 5);
%     
%     % Assign cells to clusters
%     clustered_cells = cell(numClusters, 1);
%     for k = 1:numClusters
%         clustered_cells{k} = find(idx == k);
%         fprintf('Cluster %d: %d cells\n', k, length(clustered_cells{k}));
%     end
%     
%     % Plot histogram of SNR values
%     figure;
%     histogram(snr_values, 30);
%     xlabel('SNR'); ylabel('Number of Cells');
%     title(sprintf('SNR Distribution (Clusters = %d)', numClusters));
% end
