%
dFF_all = cell2mat(combinedTable_NM.dFF);
variances = var(dFF_all, 0, 2); % Variance across time (columns)
Q1 = prctile(dFF_all, 25, 1);
Q2 = prctile(dFF_all, 75, 1);
threshold = median(variances); % You can adjust this threshold
high_noise_idx = find(variances > threshold); % Indices of high-noise curves
low_noise_idx = find(variances <= threshold); % Indices of low-noise curves

for x = 1:size(dFF_all, 1)
 %Categorize base on curve "shape"
    curve_median = median(combinedTable_NM.dFF{x,1});
    curve_mim = min(combinedTable_NM.dFF{x,1});
    curve_max = max(combinedTable_NM.dFF{x,1});
    curve_Q1 = prctile(combinedTable_NM.dFF{x,1}, 25);
    curve_Q2 = prctile(combinedTable_NM.dFF{x,1}, 75);

%     %Categorize base on curve variance across time
%     smoothed_signal = movmean(fullCraniotomy_combinedTable_NM.dFF{x,1}, 60); % 5-point moving average
%     plot(fullCraniotomy_combinedTable_NM.dFF{x,1}); hold on;
%     plot(smoothed_signal, 'r', 'LineWidth', 2); % Overlay smoothed signal in red
%     legend('Original Signal', 'Smoothed Signal');

    snr_value = snr(dFF_all(x,:));
    variance_value = var(combinedTable_NM.dFF{x,1});

    disp(['SNR: ', num2str(snr_value)]);

end

    %%

    % Convert cell array to matrix (each row = one curve, each column = time point)
dataMatrix = cell2mat(fullCraniotomy_combinedTable_NM.dFF); 

% Perform PCA
[coeff, score, ~, ~, explained] = pca(dataMatrix);

% Plot first two principal components
scatter(score(:,1), score(:,2));
xlabel('PC1'); ylabel('PC2');
title('PCA of Waveforms');

% Explained variance plot
figure;
plot(cumsum(explained), '-o');
xlabel('Number of Components'); ylabel('Cumulative Variance Explained');
title('Variance Explained by Principal Components');


numClusters = 5; % Choose number of clusters
[idx, C] = kmeans(score(:, 1:3), numClusters); % Use first 3 PCs

% Plot clustered curves
figure; hold on;
colors = lines(numClusters);
for k = 1:numClusters
    clusterIdx = find(idx == k);
    for i = 1:length(clusterIdx)
        plot(fullCraniotomy_combinedTable_NM.dFF{clusterIdx(i), 1}, 'Color', colors(k,:));
    end
end
title('Clustered Curves');
legend(arrayfun(@(x) sprintf('Cluster %d', x), 1:numClusters, 'UniformOutput', false));
