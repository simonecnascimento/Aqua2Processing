function [idx, count_k1, count_k2] = clusterWaves(dFF_all, optimalK)
    rng(42); % Set random seed for reproducibility

    % Extract features from dFF signals
    numCells = size(dFF_all, 1);
    features = zeros(numCells, 5);  % Store extracted features
    
%     for i = 1:numCells
%         signal = dFF_all(i, :);
%     
%         % Normalize the signal (using Z-score for example)
%         normalized_signal = zscore(signal);  % Or use any other normalization method
%     
%         % Feature 1: Peak amplitude (calculated after normalization)
%         features(i, 1) = max(normalized_signal) - min(normalized_signal);
%         
%         % Feature 2: Variance (calculated after normalization)
%         features(i, 2) = var(normalized_signal);
%         
%         % Feature 3: Rise time (time to peak)
%         [~, peakIdx] = max(normalized_signal);
%         features(i, 3) = peakIdx;
%         
%         % Feature 4: Decay time (time from peak to half max)
%         halfMax = max(normalized_signal) / 2;
%         decayIdx = find(normalized_signal(peakIdx:end) < halfMax, 1, 'first') + peakIdx;
%         if isempty(decayIdx)
%             decayIdx = length(normalized_signal);
%         end
%         features(i, 4) = decayIdx - peakIdx;
%         
%         % Feature 5: FFT dominant frequency
%         fft_signal = abs(fft(normalized_signal));
%         [~, maxFreqIdx] = max(fft_signal(2:end));  % Ignore DC component
%         features(i, 5) = maxFreqIdx;
%     end
    
    for i = 1:numCells
        signal = dFF_all(i, :);
        
        % Feature 1: Peak amplitude
        features(i, 1) = max(signal) - min(signal);
        
        % Feature 2: Variance
        features(i, 2) = var(signal);
        
        % Feature 3: Rise time (time to peak)
        [~, peakIdx] = max(signal);
        features(i, 3) = peakIdx;
        
        % Feature 4: Decay time (time from peak to half max)
        halfMax = max(signal) / 2;
        decayIdx = find(signal(peakIdx:end) < halfMax, 1, 'first') + peakIdx;
        if isempty(decayIdx)
            decayIdx = length(signal);
        end
        features(i, 4) = decayIdx - peakIdx;
        
        % Feature 5: FFT dominant frequency
        fft_signal = abs(fft(signal));
        [~, maxFreqIdx] = max(fft_signal(2:end)); % Ignore DC component
        features(i, 5) = maxFreqIdx;
    end
    
    % Normalize features
    features = normalize(features);
    
    % Apply PCA for dimensionality reduction
    [coeff, score] = pca(features);
    reducedFeatures = score(:, 1:2); % Use first 2 principal components

    % Perform clustering
    idx = kmeans(reducedFeatures, optimalK, 'Replicates', 5, 'MaxIter', 300, 'Start', 'plus');
    count_k1 = sum(idx == 1);
    count_k2 = sum(idx == 2);

    % --- Compute Cluster Descriptions ---
    clusterDescriptions = struct;
    clusterMedians = zeros(optimalK, 5);
    clusterIQRs = zeros(optimalK, 5);

    for c = 1:optimalK
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

    % Plot Bar Graphs for Cluster Features
    featureNames = {'Peak Amplitude', 'Variance', 'Time to Peak', 'Total Variation', 'FFT Peak Amplitude'};
    figure;
    for f = 1:5
        subplot(2,3,f);
        bar(clusterMedians(:, f)); hold on;
        errorbar(1:optimalK, clusterMedians(:, f), clusterIQRs(:, f), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        title(featureNames{f});
        xlabel('Cluster');
        ylabel(featureNames{f});
        grid on;
    end  
    sgtitle('Cluster Feature Analysis: Median ± IQR');
    
    % Plot PCA clusters
    % Define custom colors: Orange ([1 0.5 0]) and Green ([0 0.6 0])
    customColors = [1 0.5 0; 0 0.6 0]; 
    figure;
    % Scatter plot with colors based on idx
    g = gscatter(reducedFeatures(:,1), reducedFeatures(:,2), idx, customColors, 'o', 8);
    % Fill the markers
    for i = 1:length(g)
        g(i).MarkerFaceColor = g(i).Color; % Set fill color to match the edge color
    end
    % Set legend labels manually
    legend({'Cluster 1', 'Cluster 2'}, 'Location', 'best');
    %title('Clustered Cells (PCA)');
    xlabel('PC1'); ylabel('PC2');
    box off;
    
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

    %Normalize data
    data_zscore = zscore(dFF_all, 0, 2);  % Normalize each row (cell) along time axis
    zscore_min = min(data_zscore(:));  % Global min
    zscore_max = max(data_zscore(:));  % Global max

    % Plot waveforms by cluster
    figure;
    for k = 1:optimalK
        subplot(optimalK, 1, k);
        data = data_zscore(idx == k, :)';
        plot(data, 'LineWidth', 1);
        title(sprintf('Cluster %d', k));
        xlabel('Time'); ylabel('dF/F zScore');
        ylim([zscore_min, zscore_max]);  
        xlim([0,900]);
    end

    % Plot heatmap by cluster
    for k = 1:optimalK
        data = data_zscore(idx == k, :);
        data = data(1:count_k2, :); 
        % Use imagesc instead of heatmap
        figure;
        imagesc(data); % or imagesc(data) if normalization is not needed
        caxis([zscore_min, zscore_max]);  % Normalize the color scale
        colormap parula;
        colorbar;
        title(sprintf('Cluster %d Heatmap dFF zscore', k));
        xlabel('Time');
        ylabel('Cell #');
    end
end
