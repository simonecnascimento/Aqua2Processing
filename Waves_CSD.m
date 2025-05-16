%% create clusters based on preCSD x duringCSD and preCSD x postCSD response
clear all;

experiment = input('CSD or BIBN-CSD?: ', 's');  % 's' means input as string
durationCSDmin = input('Enter duration of duringCSD in minutes (1 or 2): ');

if strcmp(experiment, 'CSD')
    load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick\0.49resolution_correct\AQuA2_data_fullCraniotomy_CSD.mat');
    if durationCSDmin == 1
        toRemove = [22,23,40,80,81,90,100,129,143,144,145,214,215,216,233,248,249,250];
        combinedTable_complete(toRemove, :) = [];
    elseif durationCSDmin == 2
        toRemove = [22,23,40,80,81,100,129,144,145,214,215,216,233,248,249,250];
        combinedTable_complete(toRemove, :) = [];
    end
elseif strcmp(experiment, 'BIBN-CSD')
    if durationCSDmin == 1
        load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\BIBN\AQuA2_data_fullCraniotomy_BIBN-CSD_1min-duringCSD.mat');
        toRemove = 43;
        combinedTable_complete(toRemove, :) = [];
    elseif durationCSDmin == 2
        load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\BIBN\AQuA2_data_fullCraniotomy_BIBN-CSD_2min-duringCSD.mat');
        toRemove = 43; %same as 1min-duringCSD
        combinedTable_complete(toRemove, :) = [];
    end
end

%% Cluster 9x2 (up, none, down - acute & chronic) - CLUSTER BY EVENT RATE

% Extract cell types and cluster CSD events
cellTypes = combinedTable_complete{:, 13};
cells = ["cells"; "perivascular"; "non-perivascular"];
[eventRate_directions, eventRate_directionCounts, eventRate_directionLabels, eventRate_clusterID, rates_mHz] = clusterWaves_CSD(eventHz_byCell, cellTypes);
eventRate_clusters = [eventRate_directions; eventRate_directionCounts];
eventRate_clusters = [eventRate_clusters, cells];

combinedTable_clusters = addvars(combinedTable_complete, eventRate_clusterID, 'NewVariableNames', 'eventRate_clusterID');

% Chi-square - comparing the distribution between perivascular vs non-perivascular for each cluster
[p_values_byCluster, stats_byCluster] = chi_square_plot(eventRate_directionCounts);

% Run Chi-square test of independence - testing whether cluster and cell type are independent
[~, p_values_All, stats_All] = chi2gof_from_table(eventRate_directionCounts);

% Chi-square test by ACUTE response
[p_values_Acute, stats_Acute] = chi_square_by_acuteResponse(eventRate_directionCounts, 3);

% Chi-square test by CHRONIC response
[p_values_Chronic, stats_Chronic] = chi_square_by_chronicResponse(eventRate_directionCounts);

% Plot heatmap of conditional probabilities (3x3)
[All_matrix, P_matrix, NP_matrix] = plotConditionalProbabilities(eventRate_directionCounts, outputDir);

% Distributions 
% Plot distribution Cell type x Cluster
[clustersPercentage, roundedData] = plotClusterDistributionByCellType(eventRate_clusters, eventRate_directions, cells, outputDir);

% Plot distribution clusters x FOV
plotFOVDistributionByCluster(combinedTable_clusters, outputDir);

% Plot distribution FOV x cluster
plotClusterDistributionByFOV(combinedTable_clusters, outputDir);

%create heatmap
plotClusterHeatmap_FOVnormalized(combinedTable_clusters, [1:9], sortedFileNames, outputDir); %[1,2, 5:9] BIBN by FOV
plotClusterHeatmap_zscore(combinedTable_clusters, [2,1,4,7,3,6,9], outputDir);

% stacked bars
probs = compute_Clusterwise_probabilities(eventRate_directionCounts);
plotStackedDirectionGroups(eventRate_directionCounts, outputDir)

%% CLUSTER BY dFF

[dFF_directions, dFF_directionCounts, dFF_directionLabels, dFF_clusterID, phaseParams_dFF] = clusterParams_CSD(paramTables_allPhases, cellTypes);
dFF_clusters = [dFF_directions; dFF_directionCounts];
dFF_clusters = [dFF_clusters, cells];

combinedTable_clusters = addvars(combinedTable_clusters, dFF_clusterID, 'NewVariableNames', 'dFF_clusterID');

% Convert vectors to tables
T_eventRate = table(eventRate_clusterID, 'VariableNames', {'eventRate_clusterID'});
T_dFF = table(dFF_clusterID, 'VariableNames', {'dFF_clusterID'});

clustersTable = [T_eventRate, T_dFF, phaseParams_dFF];

% chronic increase
chronicIncrease_cluster = [];
for cellCluster = 1:size(phaseParams_dFF, 1)
    clusterID = clustersTable.eventRate_clusterID(cellCluster);
    
    if ismember(clusterID, [1 4 7])
        currentCell_dFF = phaseParams_dFF(cellCluster, :);  % row vector
        chronicIncrease_cluster = [chronicIncrease_cluster; currentCell_dFF];  % vertically concatenate
    end
end

% acute increase
acuteIncrease_cluster = [];
for cellCluster = 1:size(phaseParams_dFF, 1)
    clusterID = clustersTable.eventRate_clusterID(cellCluster);
    
    if ismember(clusterID, 2)
        currentCell_dFF = phaseParams_dFF(cellCluster, :);  % row vector
        acuteIncrease_cluster = [acuteIncrease_cluster; currentCell_dFF];  % vertically concatenate
    end
end


% 
% figure;
% for i = 1:3
%     % Extract the column, which is a cell array of numeric scalars or empty
%     cellColumn = allCells_cluster{:, i};  % this is a cell array
%     
%     % Convert the cell array to a numeric array by unwrapping each element
%     % Assume each cell contains a scalar numeric value or empty, replace empty with 0
%     numericColumn = zeros(size(cellColumn));
%     for k = 1:numel(cellColumn)
%         if isempty(cellColumn{k})
%             numericColumn(k) = 0;
%         else
%             numericColumn(k) = cellColumn{k};
%         end
%     end
%     
%     subplot(1, 3, i);
%     histogram(numericColumn);
% end
% 
% medians = zeros(1, 3);
% iqrs = zeros(1, 3);
% 
% for i = 1:3
%     cellColumn = allCells_cluster{:, i};  % cell array
%     numericColumn = zeros(size(cellColumn));
%     
%     for k = 1:numel(cellColumn)
%         if isempty(cellColumn{k})
%             numericColumn(k) = 0;
%         else
%             numericColumn(k) = cellColumn{k};
%         end
%     end
%     
%     medians(i) = median(numericColumn);
%     iqrs(i) = iqr(numericColumn);
% end
% 
% % Display results
% for i = 1:3
%     fprintf('Column %d: Median = %.4f, IQR = %.4f\n', i, medians(i), iqrs(i));
% end




%% Cluster 3x2 (up, none, down - acute & chronic)

cellTypes = combinedTable_complete{:, 13};
cells = ["cells"; "P"; "NP"];
[directions, acuteLabels, chronicLabels, acuteIDs, chronicIDs, acuteCounts, chronicCounts, rates_mHz] = clusterWaves_CSD_acute_chronic(eventHz_byCell, cellTypes);

clustersAcute = [directions; acuteCounts];
clustersAcute = [clustersAcute, cells];
clustersChronic = [directions; chronicCounts];
clustersChronic = [clustersChronic, cells];

combinedTable_clustersAcute = addvars(combinedTable_complete, acuteIDs, 'NewVariableNames', 'acuteIDs');
combinedTable_clustersChronic = addvars(combinedTable_complete, chronicIDs, 'NewVariableNames', 'chronicIDs');

% Chi-square - comparing the distribution between perivascular vs non-perivascular for each cluster
[p_values_byCluster_acute, stats_byCluster_acute] = chi_square_plot(acuteCounts);
[p_values_byCluster_chronic, stats_byCluster_chronic] = chi_square_plot(chronicCounts);

% Run Chi-square test of independence - testing whether cluster and cell type are independent
[~, p_values_All_acute, stats_All_acute] = chi2gof_from_table(acuteCounts);
[~, p_values_All_chronic, stats_All_chronic] = chi2gof_from_table(chronicCounts);

% Distribution Cell type x Cluster 
[clustersPercentage_acute, roundedData_acute] = plotClusterDistributionByCellType(clustersAcute, directions, cells, outputDir);
figure; bar(roundedData_acute,'stacked','DisplayName','roundedData_acute')
legend(directions); xticklabels(cells(2:3)); title('Acute response by cell type')

[clustersPercentage_chronic, roundedData_chronic] = plotClusterDistributionByCellType(clustersChronic, directions, cells, outputDir);
figure; bar(roundedData_chronic,'stacked','DisplayName','roundedData_chronic')
legend(directions); xticklabels(cells(2:3)); title('Chronic response by cell type')

% Heatmaps ACUTE vs CHRONIC
plotClusterHeatmap_zscore(combinedTable_clustersAcute, 1:2);
plotClusterHeatmap_zscore(combinedTable_clustersChronic, 1:2);
 
% Compare baseline between chronic clusters %1,4,7x3,6,9
preHz = str2double(eventHz_byCell(:, 2));
preHz_cluster = [preHz, combinedTable_clusters.clusterID];
compareBaselineRates(preHz, acuteIDs, chronicIDs)


%%
load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick\0.49resolution_correct\AQuA2_data_fullCraniotomy_CSD.mat')
cellLocation_indices = combinedTable{:,13};
redLabel_indices = combinedTable{:,14};

% By cell type
perivascular_indices = cellLocation_indices == 0; 
nonPerivascular_indices = cellLocation_indices == 2;
combinedTable_perivascular = combinedTable(perivascular_indices,:);
combinedTable_nonPerivascular = combinedTable(nonPerivascular_indices,:);

% dFF
dFF_all = cell2mat(combinedTable.dFF); % Convert to matrix
dFF_perivascular = cell2mat(combinedTable_perivascular.dFF);
dFF_nonPerivascular = cell2mat(combinedTable_nonPerivascular.dFF);

Fs = 1.03; % Sampling frequency

%%
combinedTable_sorted = sortrows(combinedTable, 'Max dFF', 'descend');
dFF_all_sorted = cell2mat(combinedTable_sorted.dFF); % Convert to matrix

dFF_all_sorted_preCSD = dFF_all_sorted(:, 1:1854);
dFF_all_sorted_duringCSD = dFF_all_sorted(:, 1855:1917); %60sec
dFF_all_sorted_postCSD = dFF_all_sorted(:, 1918:3772);
%dFF_all_sorted(38,:) = [];
%dFF_check = dFF_all_sorted(:, 1700:2000);

combinedTable_perivascular_sorted = sortrows(combinedTable_perivascular, 'Max dFF', 'descend');
dFF_perivascular_sorted = cell2mat(combinedTable_perivascular_sorted.dFF);

combinedTable_nonPerivascular_sorted = sortrows(combinedTable_nonPerivascular, 'Max dFF', 'descend');
dFF_nonPerivascular_sorted = cell2mat(combinedTable_nonPerivascular_sorted.dFF);

sortedMatrix = sortrows(dFF_all_sorted, 1855, "descend");

%% plot heatmap for quick visualization

minData = min(sortedMatrix(:)); 
maxData = max(sortedMatrix(:));  

% all
figure;
imagesc(sortedMatrix); 
caxis([minData, maxData]); 
colormap(flipud(gray));
colorbar;
% sorte
figure;
imagesc(dFF_perivascular_sorted);
caxis([minData, maxData]); 
colormap(flipud(gray));
colorbar;

% perivascular 
figure;
imagesc(dFF_perivascular); 
caxis([minData, maxData]); 
colormap(flipud(gray));
colorbar;
% sorte
figure;
imagesc(dFF_perivascular_sorted);
caxis([minData, maxData]); 
colormap(flipud(gray));
colorbar;

% non-perivascular 
figure;
imagesc(dFF_nonPerivascular);
caxis([minData, maxData]); 
colormap(flipud(gray));
colorbar;
% sorted
figure;
imagesc(dFF_nonPerivascular_sorted(1:65,:)); 
caxis([minData, maxData]); 
colormap(flipud(gray));
colorbar;


%% Plot selected waves
plotSelectedWaves_together(combinedTable, [88, 265, 330, 451], 'my_figure.eps');

%% Categorize waves based on frequency, duration and amplitude of events
cellTable = categorize_waves(combinedTable);

%% Plot all curves
%saveFolder = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\dFF_waves';
plot_dFF_waves(combinedTable, saveFolder, false);

%% Visualize FFT - Fast Fourier Transformation

% for each cell
for w = 1:2 %size(dFF_all, 1)
    plotFFT(dFF_all(w, :), Fs);
end

for w = 1:5 %size(dFF_all, 1)
    normalized_signal = zscore(dFF_all(w,:));
    plotFFT(normalized_signal, Fs);
end

% for all cells
plotFFT_all(dFF_all, Fs);  % dFF_all is your matrix of signals for all cells

%% 1.Cluster based on wave shape

maxClusters = 10; % Test up to 10 clusters
optimalK = findOptimalClusters(dFF_all, maxClusters);
[idx, count_k1, count_k2]  = clusterWaves(dFF_all, optimalK);

%% 2. Cluster based on SNR

% Compute SNR
snr_ma = computeSNR(dFF_all_sorted, 'ma', 5);  % Using a window size of 5
snr_fft = computeSNR(dFF_all_sorted, 'fft', 0.2);  % Using 20% of frequencies as the cutoff

maxClusters = 10; % Test up to 10 clusters
[idx, clustered_cells, optimalK, count_k1, count_k2] = clusterSNR(snr_ma, maxClusters, dFF_all_sorted);  % Use snr_ma or snr_fft

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
[optimal_k, idx, features, cluster_means, count_k1, count_k2] = clusterCa2Signals(dFF_all_sorted, Fs, poly_order, frame_size, maxClusters);



%% check differences between clusters
darkPurple = [0.5, 0, 0.5];  % Dark purple
darkGreen = [0, 0.5, 0];  
% SNR differences
snr_1 = computeSNR(dFF_all_sorted(idx == 1,:), 'ma', 5); %'ma', 5 OR 'fft', 0.2
snr_2 = computeSNR(dFF_all_sorted(idx == 2,:), 'ma', 5);
% Median and IQR
med1 = median(snr_1); 
med2 = median(snr_2);
iqr1 = iqr(snr_1); 
iqr2 = iqr(snr_2);
% Mann-Whitney U test
[p, h, stats] = ranksum(snr_1, snr_2);
fprintf('Mann-Whitney p = %.4f\n', p);
% Plot
figure;
bar(1, med1, 'FaceColor', darkPurple); hold on;
bar(2, med2, 'FaceColor', darkGreen);
%bar([med1, med2], 'FaceColor', [0.2 0.6 0.8]); %'Color', colors
% Set x-axis labels
set(gca, 'XTick', [1 2], 'XTickLabel', {'Cluster 1', 'Cluster 2'});
ylabel('Median SNR');
% Error bars (IQR)
hold on;
errorbar([1, 2], [med1, med2], [iqr1, iqr2], '.k', 'LineWidth', 1.5);
% Annotate p-value
y_max = max([med1 + iqr1, med2 + iqr2]);
text(1.5, y_max * 1.05, sprintf('Mann-Whitney p = %.4f', p),'HorizontalAlignment', 'center', 'FontSize', 12);
hold off; box off;

% cell location
wave_location = zeros(size(idx,1), 2); % Preallocate for efficiency
for cellIdx = 1:size(idx,1)
    wave_location(cellIdx, :) = [idx(cellIdx), combinedTable_sorted{cellIdx, "Cell location (0,perivascular;1,adjacent;2,none)"}];
end
tbl = chiSquaredTestForAssociation(wave_location);

% Calculate the column sums
colSums = sum(tbl);
% Convert each value to a percentage of its respective column
percentageDataColumn = (tbl ./ colSums) * 100;

% Calculate the row sums
rowSums = sum(tbl,2);
% Convert each value to a percentage of its respective column
percentageDataRow = (tbl ./ rowSums) * 100;

% other features
%combinedTable_NM(105,:) = [];
clusterFeatures = table(idx, ...
    combinedTable_sorted.("Area(um2)"), ...
    combinedTable_sorted.Perimeter, ...
    combinedTable_sorted.Circularity, ...
    combinedTable_sorted.("Max dFF"), ...
    combinedTable_sorted.("dFF AUC"),...
    combinedTable_sorted.("Duration 10% to 10%"), ...
    combinedTable_sorted.("Number of Events"), ...
    'VariableNames', {'Cluster', 'Area_um2', 'Perimeter', 'Circularity', 'Max_dFF', 'dFF AUC', 'Duration_10to10', 'Num_Events'});

% Split data into two tables based on cluster assignment
clusterFeatures_cluster1 = clusterFeatures(clusterFeatures.Cluster == 1, :);
clusterFeatures_cluster2 = clusterFeatures(clusterFeatures.Cluster == 2, :);

compareClusterFeatures(clusterFeatures, clusterFeatures_cluster1, clusterFeatures_cluster2);

%number of events 
numberOfEvents_Cluster1 = clusterFeatures_cluster1(:,"Num_Events");
totalEvents_Cluster1 = sum(numberOfEvents_Cluster1{:,:});
%number of events in Hz
numberOfEvents_Cluster1_Hz = table2cell(numberOfEvents_Cluster1);
numberOfEvents_Cluster1_Hz = cell2mat(numberOfEvents_Cluster1_Hz);
numberOfEvents_Cluster1_Hz = numberOfEvents_Cluster1_Hz / 900;
numberOfEvents_Cluster1_Hz_new = numberOfEvents_Cluster1_Hz * 1000;

%number of events 
numberOfEvents_Cluster2 = clusterFeatures_cluster2(:,"Num_Events");
totalEvents_Cluster2 = sum(numberOfEvents_Cluster2{:,:});
%number of events in Hz
numberOfEvents_Cluster2_Hz = table2cell(numberOfEvents_Cluster2);
numberOfEvents_Cluster2_Hz = cell2mat(numberOfEvents_Cluster2_Hz);
numberOfEvents_Cluster2_Hz = numberOfEvents_Cluster2_Hz / 900;
numberOfEvents_Cluster2_Hz_new = numberOfEvents_Cluster2_Hz * 1000;

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

