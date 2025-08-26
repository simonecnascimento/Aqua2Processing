%% create clusters based on preCSD x duringCSD and preCSD x postCSD response
clear all;

experiment = input('CSD or BIBN-CSD?: ', 's');  % 's' means input as string
durationCSDmin = input('Enter duration of duringCSD in minutes (1 or 2): ');

if strcmp(experiment, 'CSD')
    load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick\0.49resolution_correct\AQuA2_data_fullCraniotomy_CSD.mat');
    if durationCSDmin == 1 %due to mismatch of combinedTable and eventHz_byCell, which filters out events that happen after 61min
        toRemove = [22,23,40,80,81,90,100,129,143,144,145,214,215,216,233,248,249,250];
        combinedTable_complete(toRemove, :) = [];
    elseif durationCSDmin == 2 %due to mismatch of combinedTable and eventHz_byCell, which filters out events that happen after 62min
        toRemove = [22,23,40,80,81,100,129,144,145,214,215,216,233,248,249,250];
        combinedTable_complete(toRemove, :) = [];
    end
elseif strcmp(experiment, 'BIBN-CSD')
    if durationCSDmin == 1
        load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\BIBN\AQuA2_data_fullCraniotomy_BIBN_CSD-1min-duringCSD.mat');
        toRemove = 43;
        combinedTable_complete(toRemove, :) = [];
    elseif durationCSDmin == 2
        load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\BIBN\AQuA2_data_fullCraniotomy_BIBN-CSD_2min-duringCSD.mat');
        toRemove = 43; %same as 1min-duringCSD
        combinedTable_complete(toRemove, :) = [];
    end
end

%% CLUSTER BY EVENT RATE (9)
% 
% % Extract cell types and cluster CSD events
% cellTypes = combinedTable_complete{:, 13};
% cells = ["cells"; "perivascular"; "non-perivascular"];
% 
% % 9 clusters
% [eventRate_directions, eventRate_directionCounts, eventRate_directionLabels, eventRate_clusterID, rates_mHz, rates_mHz_clusterID] = clusterWaves_eventRate9_CSD(eventHz_byCell, cellTypes);
% eventRate_clusters = [eventRate_directions; eventRate_directionCounts];
% eventRate_clusters = [eventRate_clusters, cells];
% 
% combinedTable_clusters = addvars(combinedTable_complete, eventRate_clusterID, 'NewVariableNames', 'eventRate_clusterID');
% 
% % Assuming rates_mHz is Nx3 (columns: preCSD, duringCSD, postCSD)
% combinedTable_clusters = addvars(combinedTable_clusters, ...
%     rates_mHz(:,1), rates_mHz(:,2), rates_mHz(:,3), ...
%     'NewVariableNames', {'eventRate_preCSD', 'eventRate_duringCSD', 'eventRate_postCSD'});
% 
% rates_mHz_clusterID_sorted = sortrows(rates_mHz_clusterID, 4); % descending sort by column 4
% 
% % Chi-square - comparing the distribution between perivascular vs non-perivascular for each cluster
% [p_values_byCluster, stats_byCluster] = chi_square_plot(eventRate_directionCounts);
% 
% % Run Chi-square test of independence - testing whether cluster and cell type are independent
% [~, p_values_All, stats_All] = chi2gof_from_table(eventRate_directionCounts);
% 
% % Chi-square test by ACUTE response
% [p_values_Acute, stats_Acute] = chi_square_by_acuteResponse(eventRate_directionCounts, 3);
% 
% % Chi-square test by CHRONIC response
% [p_values_Chronic, stats_Chronic] = chi_square_by_chronicResponse(eventRate_directionCounts);
% 
% % Plot heatmap of conditional probabilities (3x3)
% [All_matrix, P_matrix, NP_matrix] = plotConditionalProbabilities(eventRate_directionCounts, outputDir);
% 
% Distributions 
% Plot distribution Cell type x Cluster
[clustersPercentage, roundedData] = plotClusterDistributionByCellType(eventRate_clusters, eventRate_directions, cells, outputDir);

% Plot distribution clusters x FOV
plotFOVDistributionByCluster(combinedTable_clusters, outputDir);

% Plot distribution FOV x cluster
plotClusterDistributionByFOV(combinedTable_clusters, outputDir);
% 
% %create heatmap
% plotClusterHeatmap_FOVnormalized(combinedTable_clusters, [1:9], sortedFileNames, outputDir); %[1,2, 5:9] BIBN by FOV
% plotClusterHeatmap_zscore(combinedTable_clusters, [2,4,9], outputDir); %[1,2,3,1,4,7,3,6,9]
% 
% % stacked bars
% probs = compute_Clusterwise_probabilities(eventRate_directionCounts);
% plotStackedDirectionGroups(eventRate_directionCounts, outputDir)

%% CLUSTER by EVENT RATE (6)

% Extract cell types and cluster CSD events
cellTypes = combinedTable_complete{:, 13};
cells = ["cells"; "perivascular"; "non-perivascular"];

% 6 cluster
[eventRate_directions, eventRate_directionCounts, eventRate_directionLabels, eventRate_clusterID, rates_mHz, rates_mHz_clusterID] = clusterWaves_eventRate6_CSD(eventHz_byCell, cellTypes);
eventRate_clusters = [eventRate_directions; eventRate_directionCounts];
eventRate_clusters = [eventRate_clusters, cells];

combinedTable_clusters = addvars(combinedTable_complete, eventRate_clusterID, 'NewVariableNames', 'eventRate_clusterID');

% % Example 1: Acute Increase
% AcuteIncrease = [7 7; 3 4];
% plotChi2Residuals(AcuteIncrease, {'P','NP'}, {'Persistent Increase','Persistent Decrease'});
% 
% % Example 2: Persistent Increase
% PersistentIncrease = [7 11; 3 34];
% plotChi2Residuals(PersistentIncrease, {'P','NP'}, {'Acute Increase','Acute No Change'});
% 
% % Example 3: Persistent Decrease
% PersistentDecrease = [7 25; 4 110];
% plotChi2Residuals(PersistentDecrease, {'P','NP'}, {'Acute Increase','Acute No Change'});
% 

% max dFF during CSD
duringCSD_maxDFF = table2cell(paramTables_allPhases.Max_dFF(:,3));

% Assuming rates_mHz is Nx3 (columns: preCSD, duringCSD, postCSD)
combinedTable_clusters = addvars(combinedTable_clusters, ...
    rates_mHz(:,1), rates_mHz(:,2), rates_mHz(:,3), duringCSD_maxDFF,...
    'NewVariableNames', {'eventRate_preCSD', 'eventRate_duringCSD', 'eventRate_postCSD', 'duringCSD_maxDFF'});

duringCSD_maxDFF_all = [];
for x = 1:size(combinedTable_clusters,1)
    if combinedTable_clusters.eventRate_clusterID(x) == 2
        duringCSD_maxDFF_x = combinedTable_clusters.duringCSD_maxDFF(x);
        duringCSD_maxDFF_all = [duringCSD_maxDFF_all; duringCSD_maxDFF_x];
    end
end



% % without separating by cell type
% rate_cluster_acuteIncrease = []; %2
% rate_cluster_chronicIncrease = []; %4
% rate_cluster_chronicDecrease = []; %6
% 
% rate_cluster_AllacuteIncrease = []; %1,2,3
% rate_cluster_AllchronicIncrease = []; %1,4
% rate_cluster_AllchronicDecrease = []; %3,6
% 
% for cellCluster = 1:size(rates_mHz, 1)
%     clusterID = combinedTable_clusters.eventRate_clusterID(cellCluster);
%     currentCell_rate = rates_mHz(cellCluster, :);  % row vector
%     cellType = cellTypes(cellCluster);
% 
%     % Individual clusters
%     if clusterID == 2
%         rate_cluster_acuteIncrease = [rate_cluster_acuteIncrease; [currentCell_rate, cellType]];
%     end
%     if clusterID == 4
%         rate_cluster_chronicIncrease = [rate_cluster_chronicIncrease; [currentCell_rate, cellType]];
%     end
%     if clusterID == 6
%         rate_cluster_chronicDecrease = [rate_cluster_chronicDecrease; [currentCell_rate, cellType]];
%     end
% 
%     % Combined clusters
%     if ismember(clusterID, [1,2,3])
%         rate_cluster_AllacuteIncrease = [rate_cluster_AllacuteIncrease; [currentCell_rate, cellType]];
%     end
%     if ismember(clusterID, [1, 4])
%         rate_cluster_AllchronicIncrease = [rate_cluster_AllchronicIncrease; [currentCell_rate, cellType]];
%     end
%     if ismember(clusterID, [3, 6])
%         rate_cluster_AllchronicDecrease = [rate_cluster_AllchronicDecrease; [currentCell_rate, cellType]];
%     end
% end

% === Initialize ===
rate_cluster_acuteIncrease_peri       = [];
rate_cluster_acuteIncrease_nonperi    = [];
rate_cluster_chronicIncrease_peri     = [];
rate_cluster_chronicIncrease_nonperi  = [];
rate_cluster_chronicDecrease_peri     = [];
rate_cluster_chronicDecrease_nonperi  = [];

% rate_cluster_AllacuteIncrease_peri    = [];
% rate_cluster_AllacuteIncrease_nonperi = [];
% rate_cluster_AllchronicIncrease_peri  = [];
% rate_cluster_AllchronicIncrease_nonperi = [];
% rate_cluster_AllchronicDecrease_peri  = [];
% rate_cluster_AllchronicDecrease_nonperi = [];

% === Loop through each cell ===
for cellCluster = 1:size(rates_mHz, 1)
    clusterID        = combinedTable_clusters.eventRate_clusterID(cellCluster);
    currentCell_rate = rates_mHz(cellCluster, :);  % row vector
    cellType         = cellTypes(cellCluster);     % 0 = peri, 2 = non-peri

    % --- Individual clusters ---
    if clusterID == 2
        if cellType == 0
            rate_cluster_acuteIncrease_peri = [rate_cluster_acuteIncrease_peri; currentCell_rate];
        elseif cellType == 2
            rate_cluster_acuteIncrease_nonperi = [rate_cluster_acuteIncrease_nonperi; currentCell_rate];
        end
    elseif clusterID == 4
        if cellType == 0
            rate_cluster_chronicIncrease_peri = [rate_cluster_chronicIncrease_peri; currentCell_rate];
        elseif cellType == 2
            rate_cluster_chronicIncrease_nonperi = [rate_cluster_chronicIncrease_nonperi; currentCell_rate];
        end
    elseif clusterID == 6
        if cellType == 0
            rate_cluster_chronicDecrease_peri = [rate_cluster_chronicDecrease_peri; currentCell_rate];
        elseif cellType == 2
            rate_cluster_chronicDecrease_nonperi = [rate_cluster_chronicDecrease_nonperi; currentCell_rate];
        end
    end

%     % --- Combined clusters ---
%     if ismember(clusterID, [1, 2, 3]) % All acute increase
%         if cellType == 0
%             rate_cluster_AllacuteIncrease_peri = [rate_cluster_AllacuteIncrease_peri; currentCell_rate];
%         elseif cellType == 2
%             rate_cluster_AllacuteIncrease_nonperi = [rate_cluster_AllacuteIncrease_nonperi; currentCell_rate];
%         end
%     end
% 
%     if ismember(clusterID, [1, 4]) % All chronic increase
%         if cellType == 0
%             rate_cluster_AllchronicIncrease_peri = [rate_cluster_AllchronicIncrease_peri; currentCell_rate];
%         elseif cellType == 2
%             rate_cluster_AllchronicIncrease_nonperi = [rate_cluster_AllchronicIncrease_nonperi; currentCell_rate];
%         end
%     end
% 
%     if ismember(clusterID, [3, 6]) % All chronic decrease
%         if cellType == 0
%             rate_cluster_AllchronicDecrease_peri = [rate_cluster_AllchronicDecrease_peri; currentCell_rate];
%         elseif cellType == 2
%             rate_cluster_AllchronicDecrease_nonperi = [rate_cluster_AllchronicDecrease_nonperi; currentCell_rate];
%         end
%     end
end

% Initialize group labels
eventRateGroup = strings(height(combinedTable_clusters), 1);
% Assign group based on event rate
eventRateGroup(combinedTable_clusters.eventRate_preCSD == 0) = "0";
eventRateGroup(combinedTable_clusters.eventRate_preCSD > 0 & combinedTable_clusters.eventRate_preCSD <= 1) = "0-1";
eventRateGroup(combinedTable_clusters.eventRate_preCSD > 1) = ">1";
% Convert to categorical
combinedTable_clusters.eventRateGroup = categorical(eventRateGroup, {'0', '0-1', '>1'});
% Use groupcounts on the new categorical variable
groupCounts = groupcounts(combinedTable_clusters.eventRateGroup);
total = sum(groupCounts);
groupPercents = 100 * groupCounts / total;

% Plot heatmap of conditional probabilities (3x3)
[All_matrix, P_matrix, NP_matrix] = plotCluster6Distribution(eventRate_directionCounts, outputDir);

%% CLUSTER by EVENT COUNT (6)

% cellTypes = combinedTable_complete{:, 13};
% cells = ["cells"; "perivascular"; "non-perivascular"];
% [eventCounts_directions, eventCounts_directionCounts, eventCounts_directionLabels, counts_clusterID, eventCounts, eventCounts_clusterID] = clusterWaves_eventCount6_CSD(events_byCell_all, cellTypes);
% eventCounts_clusters = [eventCounts_directions; eventCounts_directionCounts];
% eventCounts_clusters = [eventCounts_clusters, cells];
% 
% combinedTable_clusters = addvars(combinedTable_complete, counts_clusterID, 'NewVariableNames', 'counts_clusterID');
% 
% % Plot heatmap of conditional probabilities (3x3)
% [All_matrix, P_matrix, NP_matrix] = plotCluster6Distribution(eventCounts_directionCounts, outputDir);

%% CLUSTER BY dFF (9)

% cellTypes = combinedTable_complete{:, 13};
% cells = ["cells"; "perivascular"; "non-perivascular"];
% [dFF_directions, dFF_directionCounts, dFF_directionLabels, dFF_clusterID, phaseParams_dFF] = clusterParams_CSD(paramTables_allPhases, cellTypes);
% dFF_clusters = [dFF_directions; dFF_directionCounts];
% dFF_clusters = [dFF_clusters, cells];
% 
% combinedTable_clusters = addvars(combinedTable_clusters, dFF_clusterID, 'NewVariableNames', 'dFF_clusterID');
% 
% % Convert vectors to tables
% T_eventRate = table(eventRate_clusterID, 'VariableNames', {'eventRate_clusterID'});
% T_dFF = table(dFF_clusterID, 'VariableNames', {'dFF_clusterID'});
% 
% dFF_clustersTable = [T_eventRate, T_dFF, phaseParams_dFF];
% 
% 
% dFFcluster_acuteIncrease = []; %2
% dFFcluster_chronicIncrease = []; %4
% dFFcluster_chronicDecrease = []; %6
% 
% dFFcluster_AllchronicIncrease = []; %1,4
% dFFcluster_AllchronicDecrease = []; %3,6
% 
% 
% for cellCluster = 1:size(phaseParams_dFF, 1)
%     clusterID = dFF_clustersTable.eventRate_clusterID(cellCluster);
%     currentCell_dFF = phaseParams_dFF(cellCluster, :);  % row vector
%     
%     % Individual clusters
%     if clusterID == 2
%         dFFcluster_acuteIncrease = [dFFcluster_acuteIncrease; currentCell_dFF];
%     end
%     if clusterID == 4
%         dFFcluster_chronicIncrease = [dFFcluster_chronicIncrease; currentCell_dFF];
%     end
%     if clusterID == 6
%         dFFcluster_chronicDecrease = [dFFcluster_chronicDecrease; currentCell_dFF];
%     end
% 
%     % Combined clusters
%     if ismember(clusterID, [1, 4])
%         dFFcluster_AllchronicIncrease = [dFFcluster_AllchronicIncrease; currentCell_dFF];
%     end
%     if ismember(clusterID, [3, 6])
%         dFFcluster_AllchronicDecrease = [dFFcluster_AllchronicDecrease; currentCell_dFF];
%     end
% end


%% Cluster 3x2 (up, none, down - acute & chronic)

% cellTypes = combinedTable_complete{:, 13};
% cells = ["cells"; "P"; "NP"];
% [directions, acuteLabels, chronicLabels, acuteIDs, chronicIDs, acuteCounts, chronicCounts, rates_mHz] = clusterWaves_CSD_acute_chronic(eventHz_byCell, cellTypes);
% 
% clustersAcute = [directions; acuteCounts];
% clustersAcute = [clustersAcute, cells];
% clustersChronic = [directions; chronicCounts];
% clustersChronic = [clustersChronic, cells];
% 
% combinedTable_clustersAcute = addvars(combinedTable_complete, acuteIDs, 'NewVariableNames', 'acuteIDs');
% combinedTable_clustersChronic = addvars(combinedTable_complete, chronicIDs, 'NewVariableNames', 'chronicIDs');
% 
% % Chi-square - comparing the distribution between perivascular vs non-perivascular for each cluster
% [p_values_byCluster_acute, stats_byCluster_acute] = chi_square_plot(acuteCounts);
% [p_values_byCluster_chronic, stats_byCluster_chronic] = chi_square_plot(chronicCounts);
% 
% % Run Chi-square test of independence - testing whether cluster and cell type are independent
% [~, p_values_All_acute, stats_All_acute] = chi2gof_from_table(acuteCounts);
% [~, p_values_All_chronic, stats_All_chronic] = chi2gof_from_table(chronicCounts);
% 
% % Distribution Cell type x Cluster 
% [clustersPercentage_acute, roundedData_acute] = plotClusterDistributionByCellType(clustersAcute, directions, cells, outputDir);
% figure; bar(roundedData_acute,'stacked','DisplayName','roundedData_acute')
% legend(directions); xticklabels(cells(2:3)); title('Acute response by cell type')
% 
% [clustersPercentage_chronic, roundedData_chronic] = plotClusterDistributionByCellType(clustersChronic, directions, cells, outputDir);
% figure; bar(roundedData_chronic,'stacked','DisplayName','roundedData_chronic')
% legend(directions); xticklabels(cells(2:3)); title('Chronic response by cell type')
% 
% % Heatmaps ACUTE vs CHRONIC
% plotClusterHeatmap_zscore(combinedTable_clustersAcute, 1:2);
% plotClusterHeatmap_zscore(combinedTable_clustersChronic, 1:2);
%  
% % Compare baseline between chronic clusters %1,4,7x3,6,9
% preHz = str2double(eventHz_byCell(:, 2));
% preHz_cluster = [preHz, combinedTable_clusters.eventRate_clusterID];
% compareBaselineRates(preHz, acuteIDs, chronicIDs)

%% EVENT raster plot

% Process startingFrames for raster plot
allEmpty = all(cellfun(@isempty, startingFrames_byCell_all(:, 2:5)), 2);
startingFrames_byCell_all(allEmpty, :) = []; %% Remove rows where all event phases are empty

% Step 1: Flatten nested cells inside columns 2 to 5
flattenedData = startingFrames_byCell_all;

for i = 1:size(flattenedData, 1)
    for j = 2:5
        val = flattenedData{i, j};
        if iscell(val)
            % If it's a cell, convert to numeric vector
            try
                flattenedData{i, j} = cell2mat(val);
            catch
                flattenedData{i, j} = [];
            end
        end
    end
end

% Step 2: Convert to table
startingFrames_table = cell2table(flattenedData, ...
    'VariableNames', {'cellID', 'preCSD', 'duringCSD', 'postCSD', 'baseline_preCSD'});

% Convert clusterID column to a table
clusterID_table = table(combinedTable_clusters.eventRate_clusterID(:), 'VariableNames', {'eventRate_clusterID'});
dFF_table = table(combinedTable_clusters.dFF(:), 'VariableNames', {'dFF'});

% Horizontally concatenate the tables
rasterTable = [startingFrames_table, clusterID_table, dFF_table];
rasterTable_sorted = plotRasterByCluster(rasterTable, [2,4,6], experiment);

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

