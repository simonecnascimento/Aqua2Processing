function optimalK = findOptimalClusters(dFF_all, maxClusters)
    % Extract features for clustering
    numCells = size(dFF_all, 1);
    features = zeros(numCells, 5);  
    
    for i = 1:numCells
        signal = dFF_all(i, :);
        features(i, 1) = max(signal) - min(signal); % Peak amplitude
        features(i, 2) = var(signal);              % Variance
        features(i, 3) = find(signal == max(signal), 1); % Time to peak
        features(i, 4) = sum(abs(diff(signal)));   % Total signal change
        fft_signal = abs(fft(signal));             % FFT feature
        features(i, 5) = max(fft_signal(2:end));
    end

    % Normalize features
    features = normalize(features);

    % Use PCA to reduce dimensions
    [~, score] = pca(features);
    reducedFeatures = score(:, 1:2); % Use first 2 principal components

    % Compute sum of squared errors (SSE) for different cluster numbers
    sum_of_squares = zeros(maxClusters, 1);
    for k = 1:maxClusters
        [~, ~, sumd] = kmeans(reducedFeatures, k, 'Replicates', 5);
        sum_of_squares(k) = sum(sumd);
    end

%     % Plot Elbow Curve
%     figure;
%     plot(1:maxClusters, sum_of_squares, '-o', 'LineWidth', 2);
%     xlabel('Number of Clusters');
%     ylabel('Sum of Squared Distances (SSE)');
%     title('Elbow Method for Optimal Clusters');
%     grid on;

    % Find the "elbow" point using the second derivative method
    diffs = diff(sum_of_squares);
    second_diffs = diff(diffs);
    [~, optimalK] = max(abs(second_diffs));
    optimalK = optimalK + 1; % Adjust index because diff reduces size

    fprintf('Optimal number of clusters: %d\n', optimalK);

    figure;
    plot(1:maxClusters, sum_of_squares, '-o', 'LineWidth', 2);
    hold on;
    plot(optimalK, sum_of_squares(optimalK), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Number of Clusters'); ylabel('Sum of Squared Distances');
    title(sprintf('Elbow Method: Optimal Clusters = %d', optimalK));
    grid on;
    hold off;
end
