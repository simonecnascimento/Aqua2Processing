function [optimal_k, idx, features, cluster_means, count_k1, count_k2] = clusterCa2Signals(dFF_all, Fs, poly_order, frame_size, maxClusters)
    % CLUSTERCa2SIGNALS Clusters Ca2+ signals based on frequency domain features
    %
    % Inputs:
    %   dFF_all      - Matrix where each row represents a Ca2+ signal for each cell
    %   Fs           - Sampling frequency in Hz
    %   poly_orders  - Array of polynomial orders for Savitzky-Golay filter
    %   frame_sizes  - Array of frame sizes for Savitzky-Golay filter
    %   maxClusters  - Maximum number of clusters to consider for k-means
    %
    % Outputs:
    %   optimal_k    - Optimal number of clusters based on Elbow and Silhouette methods
    %   idx          - Cluster indices for each signal
    %   features     - Extracted features (dominant frequency, peak count, etc.)
    %   cluster_means - Mean signal of each cluster

    rng(42); % Set random seed for reproducibility

    % Number of signals (cells)
    numSignals = size(dFF_all, 1);

    % Pre-allocate feature matrix for clustering
    features = [];

    psd_all = [];

    % Loop through each signal (each row of dFF_all)
    for w = 1:numSignals

        % 1. Detrend the signal (remove systematic drifts like tissue drift, bleaching)
        detrended_signal = detrend(dFF_all(w, :));
        %figure; plot(dFF_all(w, :));
        %figure; plot(detrended_signal);

        % 2. Apply Savitzky-Golay filter for smoothing (window size 5, polynomial order 3)
        smoothed_signal = sgolayfilt(detrended_signal, poly_order, frame_size);
        %figure; plot(smoothed_signal);

        % 3. Fourier Transformation: Calculate the power spectrum density (PSD)
        L = length(smoothed_signal);  % Length of the signal

        % Compute FFT
        Y = fft(smoothed_signal);     
        P2 = abs(Y/L);                
        P1 = P2(1:L/2+1);            
        P1(2:end-1) = 2*P1(2:end-1);  % Adjust for single-sided spectrum

        % Frequency axis
        f = Fs * (0:(L/2)) / L;

        % Compute Power Spectrum Density (PSD)
        psd = P1.^2;
        psd_all = [psd_all; psd];

        % 4. Peak Detection (using MPP and HPA ruleset)
        MPP = 0.1;  % For example, 10% of the max amplitude as MPP
        HPA = max(smoothed_signal) - min(smoothed_signal);

        smoothed_signal = smoothed_signal(:);
        % Remove NaN and Inf values
        smoothed_signal = smoothed_signal(~isnan(smoothed_signal) & ~isinf(smoothed_signal));

        % Find peaks with a minimum prominence
        [peaks, locs] = findpeaks(smoothed_signal, 'MinPeakProminence', MPP * HPA);

        % 5. Frequency Domain Features (using dominant frequency)
        [~, idx] = max(psd); 
        dominant_frequency = f(idx);

        % 6. Extract features for clustering (e.g., dominant frequency, peak count, etc.)
        features = [features; dominant_frequency, length(peaks)];
    end

    % 7. Determine optimal number of clusters using Elbow Method and Silhouette Analysis
    % Calculate inertia (sum of squared distances to the nearest cluster center) for different k
    inertia = zeros(maxClusters, 1);
    silhouette_scores = zeros(maxClusters, 1);

    for k = 1:maxClusters
        % Apply k-means clustering
        [idx, C, sumd] = kmeans(features, k);

        % Compute inertia (within-cluster sum of squared distances)
        inertia(k) = sum(sumd);

        % Compute silhouette score
        silhouette_scores(k) = mean(silhouette(features, idx));
    end

    % Elbow Method: Plot Inertia vs Number of Clusters
    figure;
    subplot(1, 2, 1);
    plot(1:maxClusters, inertia, '-o', 'LineWidth', 2);
    title('Elbow Method');
    xlabel('Number of Clusters');
    ylabel('Inertia');
    grid on;

%     % Silhouette Method: Plot Silhouette Scores vs Number of Clusters
%     subplot(1, 2, 2);
%     plot(1:maxClusters, silhouette_scores, '-o', 'LineWidth', 2);
%     title('Silhouette Method');
%     xlabel('Number of Clusters');
%     ylabel('Average Silhouette Score');
%     grid on;

    % Choose the optimal number of clusters based on the methods above
    [~, optimal_k_elbow] = min(diff(inertia)); % Elbow method (could also use other methods)
    %optimal_k_silhouette = find(silhouette_scores == max(silhouette_scores)); % Silhouette method

    % You can print or visualize the optimal k chosen from either method
    disp(['Optimal number of clusters (Elbow method): ', num2str(optimal_k_elbow + 1)]);
    %disp(['Optimal number of clusters (Silhouette method): ', num2str(optimal_k_silhouette)]);

    % Apply clustering with the optimal k (e.g., from the elbow method)
    optimal_k = optimal_k_elbow + 1; % Adjust to account for index difference

    % Apply k-means clustering with the optimal k
    [idx, C] = kmeans(features, optimal_k,  'Replicates', 5, 'MaxIter', 300, 'Start', 'plus');
    count_k1 = sum(idx == 1);
    count_k2 = sum(idx == 2);

    % Define custom colors
    darkPurple = [0.5, 0, 0.5];  % Dark purple
    darkGreen = [0, 0.5, 0];     % Dark green
    colors = [darkPurple; darkGreen];  % Colors for Cluster 1 and 2

    % Visualize Clusters
    figure;
    h = gscatter(features(:, 1), features(:, 2), idx, colors, '.', 20);
    xlabel('Dominant Frequency (Hz)');
    ylabel('Number of Peaks');
    legend(h, {'Cluster 1', 'Cluster 2'}, 'Location', 'best'); legend boxoff;
    set(gca, 'TickDir', 'out');
    grid off;
    box off;

    fs = 1.03;  % <-- Replace with your actual sampling rate (Hz)
    numSamples = size(dFF_all, 2);
    time_in_min = (0:numSamples-1) / fs / 60;
    
    % Plot the Ca2+ signals for each cluster, superimposesd
    figure;
    hold on;
    
    plotOrder = [2, 1];  % Cluster 2 behind Cluster 1
    legendHandles = gobjects(optimal_k, 1);
    legendLabels = cell(1, optimal_k);
    
    for i = 1:length(plotOrder)
        clusterIdx = plotOrder(i);
        clusterSignalsIdx = find(idx == clusterIdx);
    
        for sigIdx = clusterSignalsIdx'
            plot(time_in_min, dFF_all(sigIdx, :), 'Color', colors(clusterIdx, :), 'LineWidth', 1);
        end
    
        clusterMeanSignal = mean(dFF_all(clusterSignalsIdx, :), 1);
        h = plot(time_in_min, clusterMeanSignal, 'Color', colors(clusterIdx, :), 'LineWidth', 1); %'LineStyle', '-'
        legendHandles(clusterIdx) = h;
        legendLabels{clusterIdx} = ['Cluster ' num2str(clusterIdx)];
    end
    
    %title('Ca2+ Signal Curves for Each Cluster');
    xlabel('Time (min)');
    ylabel('dF/F');
    ylim([-0.5,6.5]);
    legend(legendHandles, legendLabels, 'Location', 'best'); legend boxoff;
    grid off;
    hold off;

    % Output the final clustering information
    cluster_means = C;  % Cluster centers (mean of each cluster)

    % Plot waveforms by cluster
    % Find the global min and max
    minData = min(dFF_all(:));  % Global min
    maxData = max(dFF_all(:));  % Global max

    figure;
    for k = 1:optimal_k
        subplot(optimal_k, 1, k);
        plot(dFF_all(idx == k, :)', 'LineWidth', 0.5, 'Color', colors(k,:));
            title(sprintf('Cluster %d', k));
        xlabel('Time (min)'); ylabel('\DeltaF/F');
        ylim([minData, maxData]);
        xticks = linspace(1, numSamples, 6);  % Adjust number of ticks as needed
        xticklabels = round(time_in_min(round(xticks)), 0);  % Show 2 decimal places
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
        %xlim([0,900]);
        xlim([1, numSamples]);
        box off;
    end

    % Plot heatmap by cluster
    dFF_cluster1 = dFF_all(idx == 1, :);
    dFF_cluster2 = dFF_all(idx == 2, :);
    for k = 1:optimal_k
        data = dFF_all(idx == k, :);
        if count_k1 > count_k2
            data = data(1:count_k2, :); 
        else
            data = data(1:count_k1, :); 
        end
         % Use imagesc instead of heatmap
        figure;
        imagesc(data); % or imagesc(data) if normalization is not needed
        caxis([minData, maxData]); % Normalize the color scale
        colormap(flipud(gray));
        colorbar;
        xticks = linspace(1, numSamples, 6);  % Adjust number of ticks as needed
        xticklabels = round(time_in_min(round(xticks)), 0);  % Show 2 decimal places
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
        title(sprintf('Cluster %d Heatmap dFF', k));
        xlabel('Time (min)');
        ylabel('Cell #');
        %exportgraphics(gcf, sprintf('Cluster_%d_Heatmap.eps', k), 'ContentType', 'vector');
    end
end