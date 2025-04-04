function [idx, clusterDescriptions] = clusterAndDescribe(dFF_all, numClusters)
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
        features(i, 5) = max(fft_signal(2:end));   % Peak FFT amplitude
    end

    % Normalize features
    features = normalize(features);

    % Perform k-means clustering
    [idx, C] = kmeans(features, numClusters, 'Replicates', 5);

    % Store cluster descriptions
    clusterDescriptions = struct;
    
    % Initialize arrays for plotting
    clusterMedians = zeros(numClusters, 5);
    clusterIQRs = zeros(numClusters, 5);
    
    for c = 1:numClusters
        clusterCells = find(idx == c);
        clusterData = features(clusterCells, :);

        % Compute statistics
        clusterDescriptions(c).ClusterNumber = c;
        clusterDescriptions(c).NumCells = length(clusterCells);
        clusterDescriptions(c).MedianAmplitude = median(clusterData(:,1));
        clusterDescriptions(c).IQRAmplitude = iqr(clusterData(:,1));
        clusterDescriptions(c).MedianTimeToPeak = median(clusterData(:,3));
        clusterDescriptions(c).IQRTimeToPeak = iqr(clusterData(:,3));
        clusterDescriptions(c).MedianTotalVariation = median(clusterData(:,4));
        clusterDescriptions(c).IQRTotalVariation = iqr(clusterData(:,4));
        clusterDescriptions(c).MedianFFTAmplitude = median(clusterData(:,5));
        clusterDescriptions(c).IQRFFTAmplitude = iqr(clusterData(:,5));
        
        % Store for plotting
        clusterMedians(c, :) = median(clusterData);
        clusterIQRs(c, :) = iqr(clusterData);
        
        % Display cluster summary
        fprintf('Cluster %d: %d cells\n', c, length(clusterCells));
        fprintf('  - Median Amplitude: %.2f (IQR: %.2f)\n', clusterDescriptions(c).MedianAmplitude, clusterDescriptions(c).IQRAmplitude);
        fprintf('  - Median Time to Peak: %.2f (IQR: %.2f)\n', clusterDescriptions(c).MedianTimeToPeak, clusterDescriptions(c).IQRTimeToPeak);
        fprintf('  - Median Total Variation: %.2f (IQR: %.2f)\n', clusterDescriptions(c).MedianTotalVariation, clusterDescriptions(c).IQRTotalVariation);
        fprintf('  - Median FFT Peak Amplitude: %.2f (IQR: %.2f)\n\n', clusterDescriptions(c).MedianFFTAmplitude, clusterDescriptions(c).IQRFFTAmplitude);
    end

    % Plot bar graphs with error bars (Median ± IQR)
    featureNames = {'Peak Amplitude', 'Variance', 'Time to Peak', 'Total Variation', 'FFT Peak Amplitude'};
    
    figure;
    for f = 1:5
        subplot(2,3,f);
        bar(clusterMedians(:, f)); hold on;
        errorbar(1:numClusters, clusterMedians(:, f), clusterIQRs(:, f), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        title(featureNames{f});
        xlabel('Cluster');
        ylabel(featureNames{f});
        grid on;
    end
    
    sgtitle('Cluster Feature Analysis: Median ± IQR');

end
