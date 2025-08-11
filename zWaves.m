clear all;

load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\AQuA2_data_fullCraniotomy_features_baseline.mat')
load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\multinucleated cells\multinucleatedCells.mat');
combinedTable_multinucleated = addvars(combinedTable, multinucleatedCells.Var3, 'NewVariableNames', 'Multinucleated');
multinucleated_indices = combinedTable_multinucleated{:,17};
combinedTable_NM = combinedTable_multinucleated(multinucleated_indices == 0, :);
cellLocation_indices = combinedTable_NM{:,13};
redLabel_indices = combinedTable_NM{:,14};

% By cell type
perivascular_indices = cellLocation_indices == 0; 
nonPerivascular_indices = cellLocation_indices == 2;
perivascular = combinedTable_NM(perivascular_indices,:);
nonPerivascular = combinedTable_NM(nonPerivascular_indices,:);

%% Plot selected waves
plotSelectedWaves_together(combinedTable_NM, [88, 265, 330, 451], 'my_figure.eps');

%% Categorize waves based on frequency, duration and amplitude of events
cellTable = categorize_waves(combinedTable_NM);

%% Plot all curves
saveFolder = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\dFF_waves';
plot_dFF_waves(combinedTable_NM, saveFolder, false);


%% Cluster based on wave shape



%% Cluster based on SNR

%Compute SNR first (without sgolayfilt).
%Then, use FFT to analyze frequency content.
%Apply Savitzky-Golay filtering only for visualization, not for clustering.




dFF_all = cell2mat(combinedTable_NM.dFF);  % Convert to matrix
numClusters = 4;  % Adjust based on the elbow method
clusterWaves(dFF_all, numClusters);




%% Visualize FFT - Fast Fourier Transformation

Fs = 1.03; % Sampling frequency
for w = 1:size(dFF_all, 1)
    plotFFT(dFF_all(w, :), Fs);
end

%% 1. Categorize waves based on SNR

dFF_all = cell2mat(combinedTable_NM.dFF);

% Optimize Sgolay params
poly_orders = 1:6; % Test polynomial orders from 2 to 5
frame_sizes = 5:2:927; % Test odd frame sizes from 5 to 45

[best_poly_order, best_frame_size, best_snr, best_filtered_signals] = optimizeSgolayParams(dFF_all, poly_orders, frame_sizes);



poly_order = 3;
frame_range = 5:2:927;  % Must be odd values


% Find the best frame sizes
best_frame_sizes = optimizeFrameSize(dFF_all, poly_order, frame_range);



%Compute SNR
snr_values = computeSNR(dFF_all, poly_order, best_frame_sizes);

%% Cluster SNR
maxClusters = 10;  % Maximum number of clusters to test
numClusters = 4;   % Number of clusters to use
[idx, clustered_cells] = clusterSNR(snr_values, maxClusters, numClusters);

for k = 1:numClusters
    selectedCells_indices = clustered_cells{k}; % Get indices of current cluster
    fprintf('Plotting Cluster %d with %d cells...\n', k, length(selectedCells_indices));
    
    % Plot dFF traces for the current cluster
    plot_selected_dFF_waves(combinedTable_NM, sprintf('Cluster %d', k), false, selectedCells_indices);
end

%% PCA and Kmeans clustering of dFF

dFF_all = cell2mat(combinedTable_NM.dFF);

% Define parameters
numPCs = 10;        % Number of principal components
maxClusters = 10;   % Maximum number of clusters to test

% Run automatic clustering
[idx, clustered_cells, numClusters] = autoCluster_dFF_PCA(dFF_all, numPCs, maxClusters);

% Display results
fprintf('Optimal number of clusters: %d\n', numClusters);
for k = 1:numClusters
    fprintf('Cluster %d: %d cells\n', k, length(clustered_cells{k}));
end

%% Plot waves based on SNR visual

histogram(snr_values);
Q2 = prctile(snr_values, 75);

% Define threshold for high SNR
snr_threshold = 27;

% Find high SNR cells
highSNRcells = snr_values > snr_threshold;
lowSNRcells = snr_values < snr_threshold;

% Get logical indices
selectedCells_indices = find(highSNRcells); 

% Plot selected dFF traces
set(0, 'DefaultFigureWindowStyle', 'docked');
plot_selected_dFF_waves(combinedTable_NM, '', false, selectedCells_indices);

