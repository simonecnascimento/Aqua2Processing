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

dFF_all = cell2mat(combinedTable_NM.dFF); % Convert to matrix
Fs = 1.03; % Sampling frequency
%% Plot selected waves
plotSelectedWaves_together(combinedTable_NM, [88, 265, 330, 451], 'my_figure.eps');

%% Categorize waves based on frequency, duration and amplitude of events
cellTable = categorize_waves(combinedTable_NM);

%% Plot all curves
saveFolder = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\dFF_waves';
plot_dFF_waves(combinedTable_NM, saveFolder, false);

%% Visualize FFT - Fast Fourier Transformation

% for each cell
for w = 1:2 %size(dFF_all, 1)
    plotFFT(dFF_all(w, :), Fs);
end

% for all cells
plotFFT_all(dFF_all, Fs);  % dFF_all is your matrix of signals for all cells

%% 1.Cluster based on wave shape

maxClusters = 10; % Test up to 10 clusters
optimalK = findOptimalClusters(dFF_all, maxClusters);
[idx, count_k1, count_k2]  = clusterWaves(dFF_all, optimalK);

%% 2. Cluster based on SNR

dFF_all(105, :) = [];

% Compute SNR
snr_ma = computeSNR(dFF_all, 'ma', 5);  % Using a window size of 5
snr_fft = computeSNR(dFF_all, 'fft', 0.2);  % Using 20% of frequencies as the cutoff

maxClusters = 10; % Test up to 10 clusters
[idx, clustered_cells, optimalK, count_k1, count_k2] = clusterSNR(snr_ma, maxClusters, dFF_all);  % Use snr_ma or snr_fft

%% 3. Cluster based on detrend, sgolayfilt, kmeans
% https://www-science-org.ezp-prod1.hul.harvard.edu/doi/pdf/10.1126/scisignal.abe6909

% 1. Apply Savitzky-Golay filtering only for visualization, not for clustering.
% Optimize Sgolay params - takes very long time to optimize both params
poly_orders = 1:6; 
frame_sizes = 5:2:927; 
[best_poly_order, best_frame_size, best_snr, best_filtered_signals] = optimizeSgolayParams(dFF_all, poly_orders, frame_sizes);

% Check best frame size
unique_vals_frame = unique(best_frame_size);
counts = histc(best_frame_size, unique_vals_frame);
figure; bar(unique_vals_frame, counts, 'FaceColor', 'b');

% Check best poly order
unique_vals_order = unique(best_poly_order);
counts = histc(best_poly_order, unique_vals_order);
figure; bar(unique_vals_order, counts, 'FaceColor', 'b');

% 2. Cluster signals
% Define your inputs
maxClusters = 10; % Test up to 10 clusters
poly_order = 6;  % Range of polynomial orders for Savitzky-Golay filter
frame_size = 7;  % Frame sizes (odd numbers)
% Assuming 'dFF_all' is your matrix of Ca2+ signals (cells x time points)
[optimal_k, idx, features, cluster_means, count_k1, count_k2] = clusterCa2Signals(dFF_all, Fs, poly_order, frame_size, maxClusters);


%% check differences between clusters

% cell location
wave_location = zeros(size(idx,1), 2); % Preallocate for efficiency
for cellIdx = 1:size(idx,1)
    wave_location(cellIdx, :) = [idx(cellIdx), combinedTable_NM{cellIdx, "Cell location (0,perivascular;1,adjacent;2,none)"}];
end
tbl = chiSquaredTestForAssociation(wave_location);

% Calculate the column sums
colSums = sum(tbl);
% Convert each value to a percentage of its respective column
percentageData = (tbl ./ colSums) * 100;

% other features
%combinedTable_NM(105,:) = [];
clusterFeatures = table(idx, ...
    combinedTable_NM.("Area(um2)"), ...
    combinedTable_NM.Perimeter, ...
    combinedTable_NM.Circularity, ...
    combinedTable_NM.("Max dFF"), ...
    combinedTable_NM.("Duration 10% to 10%"), ...
    combinedTable_NM.("Number of Events"), ...
    'VariableNames', {'Cluster', 'Area_um2', 'Perimeter', 'Circularity', 'Max_dFF', 'Duration_10to10', 'Num_Events'});

% Split data into two tables based on cluster assignment
clusterFeatures_cluster1 = clusterFeatures(clusterFeatures.Cluster == 1, :);
clusterFeatures_cluster2 = clusterFeatures(clusterFeatures.Cluster == 2, :);

compareClusterFeatures(clusterFeatures, clusterFeatures_cluster1, clusterFeatures_cluster2);

%%
% Plot waves based on SNR (visual)
figure;
histogram(snr_ma);
Q2 = prctile(snr_ma, 75);

% Define threshold for high SNR
snr_threshold = 22;

% Find high SNR cells
highSNRcells = snr_ma > snr_threshold;
lowSNRcells = snr_ma < snr_threshold;

% Get logical indices
selectedCells_indices = find(highSNRcells); 

% Plot selected dFF traces
set(0, 'DefaultFigureWindowStyle', 'docked');
plot_selected_dFF_waves(combinedTable_NM, '', false, selectedCells_indices);

% correlation SNR and event duration
snr_duration = zeros(size(snr_ma,1), 2);
for cellIdx = 1:size(snr_ma,1)
    snr_duration(cellIdx, :) = [snr_ma(cellIdx), combinedTable_NM{cellIdx, "Duration 10% to 10%"}];
end

%% Apply Savitzky-Golay filtering only for visualization, not for clustering.

% Optimize Sgolay params - takes very long time to optimize both params
poly_orders = 1:6; 
frame_sizes = 5:2:927; 
[best_poly_order, best_frame_size, best_snr, best_filtered_signals] = optimizeSgolayParams(dFF_all, poly_orders, frame_sizes);

% Optimize only frame range
poly_order = 3;
frame_range = 5:2:927;  % Must be odd values
best_frame_sizes = optimizeFrameSize(dFF_all, poly_order, frame_range);

