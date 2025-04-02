function [idx, clusterCenters, features] = cluster_dFF(dFF_all, numClusters)
    % CLUSTER_DFF clusters dFF traces based on waveform features.
    %
    %   [idx, clusterCenters, features] = cluster_dFF(dFF_all, numClusters)
    %
    %   Inputs:
    %       dFF_all - (cells x time points) matrix of dF/F traces.
    %       numClusters - Number of clusters for K-Means.
    %
    %   Outputs:
    %       idx - Cluster assignments for each cell.
    %       clusterCenters - Cluster center positions.
    %       features - Extracted features (used for clustering).

    [numCells, numTimePoints] = size(dFF_all);
    features = zeros(numCells, 5); % Feature matrix

    % Sampling Frequency Estimation (Assume evenly spaced time points)
    Fs = numTimePoints / 900; % From your previous info: 927 time points = 900 sec

    for i = 1:numCells
        signal = dFF_all(i, :);

        % Amplitude Features
        features(i, 1) = mean(signal);   % Mean dF/F
        features(i, 2) = max(signal);    % Max dF/F
        features(i, 3) = std(signal);    % Standard deviation

        % Frequency Feature (Dominant Frequency)
        fft_signal = abs(fft(signal));
        freq_axis = linspace(0, Fs/2, numTimePoints/2); % Half-spectrum
        [~, peakIdx] = max(fft_signal(1:numTimePoints/2));
        features(i, 4) = freq_axis(peakIdx); % Dominant frequency

        % Shape Features
        features(i, 5) = skewness(signal);  % Skewness
        features(i, 6) = kurtosis(signal); % Kurtosis
    end

    % Apply PCA for dimensionality reduction
    [coeff, score, ~] = pca(features);
    reducedFeatures = score(:, 1:3); % Keep first 3 principal components

    % Perform K-Means Clustering
    [idx, clusterCenters] = kmeans(reducedFeatures, numClusters, 'Replicates', 5);

    % Plot Clusters (2D PCA space)
    figure;
    scatter(reducedFeatures(:, 1), reducedFeatures(:, 2), 50, idx, 'filled');
    xlabel('PC1'); ylabel('PC2'); title('Clustered Cells in PCA Space');
    colorbar; grid on;

    % Histogram of SNR per cluster
    figure;
    histogram(idx, numClusters);
    xlabel('Cluster Index'); ylabel('Number of Cells');
    title('Cluster Distribution');
end
