%% Create heatmap of dFF data by experiment

%Perivascular cells
for perivascular = resultsFinal.("Cell location (0,perivascular;1,adjacent;2,none)") == 0
    dFF_perivascular = resultsFinal{perivascular,"dFF"};
end

%Convert cell array to 2D matrix
dFF_perivascular = cell2mat(dFF_perivascular);

% Define the colormap
cmap = redbluecmap; % Or any other colormap you prefer
newCmap = imresize(cmap, [128, 3]); % Increase number of colors
newCmap = min(max(newCmap, 0), 1);

% Define the color limits based on your data range
maxDff = max(dFF_perivascular(:));
minDff = min(dFF_perivascular(:));
colorLimits = [minDff maxDff]; % Using the actual minimum and maximum values of the data

% Create the heatmap with the custom colormap and color limits using imagesc
imagesc(dFF_perivascular, colorLimits);

% Apply the custom colormap
colormap(newCmap);

% Add color bar
colorbar;

% Set titles and labels
title('dFF');
xlabel('Seconds');
ylabel('Cell ID');


%% 

%Non-perivascular cells
for non_perivascular = resultsFinal.("Cell location (0,perivascular;1,adjacent;2,none)") == 2
    dFF_Nperivascular = resultsFinal{non_perivascular,"dFF"};
end

%Convert cell array to 2D matrix
dFF_Nperivascular = cell2mat(dFF_Nperivascular);


% Define the colormap
cmap = redbluecmap; % Or any other colormap you prefer
newCmap = imresize(cmap, [128, 3]); % Increase number of colors
newCmap = min(max(newCmap, 0), 1);

% Define the color limits based on your data range
maxDff = max(dFF_Nperivascular(:));
minDff = min(dFF_Nperivascular(:));
colorLimits = [minDff maxDff]; % Using the actual minimum and maximum values of the data

% Create the heatmap with the custom colormap and color limits using imagesc
imagesc(dFF_Nperivascular, colorLimits);

% Apply the custom colormap
colormap(newCmap);

% Add color bar
colorbar;

% Set titles and labels
title('dFF');
xlabel('Seconds');
ylabel('Cell ID');


%% Multiple experiments

clear all;
%for fullCraniotomy baseline data - load spreadsheet and add Var3(multinucleated cells)
load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\AQuA2_data_fullCraniotomy_features_baseline.mat')
load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\multinucleated cells\multinucleatedCells.mat');
fullCraniotomy_combinedTable = addvars(combinedTable, multinucleatedCells.Var3, 'NewVariableNames', 'Multinucleated');
fullCraniotomy_multinucleated = fullCraniotomy_combinedTable{:,17};
fullCraniotomy_combinedTable_NM = fullCraniotomy_combinedTable(fullCraniotomy_multinucleated == 0, :);
fullCraniotomy_cellLocation = fullCraniotomy_combinedTable_NM{:,13};
fullCraniotomy_redLabel = fullCraniotomy_combinedTable_NM{:,14};
% By cell type
perivascular_indices = fullCraniotomy_cellLocation == 0; % Indices of cells belonging to group 0
nonPerivascular_indices = fullCraniotomy_cellLocation == 2; % Indices of cells belonging to group 2
combinedTable_perivascular = fullCraniotomy_combinedTable(perivascular_indices,:);
combinedTable_nonPerivascular = fullCraniotomy_combinedTable(nonPerivascular_indices,:);

%% Separate by type of cell and Sort by number of events
fullCraniotomy_combinedTable_perivascular_NM_sorted = sortrows(combinedTable_perivascular, 'Number of Events', 'descend');
fullCraniotomy_combinedTable_non_perivascular_NM_sorted = sortrows(combinedTable_nonPerivascular, 'Number of Events', 'descend');

% All cells
% Extract dFF data from combinedTable, converting cell array to matrix
dFF_data_all_NM = cell2mat(fullCraniotomy_combinedTable_NM.dFF);
%perivascular cells
dFF_data_Perivascular = fullCraniotomy_combinedTable_perivascular_NM_sorted.dFF;
dFF_data_Perivascular = cell2mat(dFF_data_Perivascular);
%Non-perivascular cells
dFF_data_NonPerivascular = fullCraniotomy_combinedTable_non_perivascular_NM_sorted.dFF;
dFF_data_NonPerivascular = cell2mat(dFF_data_NonPerivascular);
dFF_data_NonPerivascular_reduced = dFF_data_NonPerivascular(1:122,:); %to match number of perivascular cells
% Find global min and max values
minVal = min([dFF_data_Perivascular(:); dFF_data_NonPerivascular(:)]);
maxVal = max([dFF_data_Perivascular(:); dFF_data_NonPerivascular(:)]);
caxis([minVal, maxVal]); % Set color axis


% Z SCORE
% Apply z-score normalization to the dFF data
dFF_data_zscore = zscore(dFF_data_all_NM, [], 2);
% Apply z-score normalization to the dFF data
dFF_data_Perivascular_zscore = zscore(dFF_data_Perivascular, [], 1);
dFF_data_NonPerivascular_zscore = zscore(dFF_data_NonPerivascular, [], 1);
dFF_data_NonPerivascular_reduced_zscore = zscore(dFF_data_NonPerivascular_reduced, [], 1);
% Find global min and max values
minVal = min([dFF_data_Perivascular_zscore(:); dFF_data_NonPerivascular_reduced_zscore(:)]);
maxVal = max([dFF_data_Perivascular_zscore(:); dFF_data_NonPerivascular_reduced_zscore(:)]);
caxis([minVal, maxVal]); % Set color axis


% Combine the data from both groups
combined_data = [dFF_data_NonPerivascular; dFF_data_Perivascular];
combined_data_zscore = [dFF_data_Perivascular_zscore; dFF_data_NonPerivascular_zscore];

% Determine the colorbar limits based on the combined data
caxis_limits_all_NM = [min(combined_data(:)), max(combined_data(:))];
caxis_limits_perivascular = [min(dFF_data_Perivascular_zscore(:)), max(dFF_data_Perivascular_zscore(:))];
caxis_limits_nonPerivascular = [min(dFF_data_NonPerivascular_zscore(:)), max(dFF_data_NonPerivascular_zscore(:))];

%% Sorted by cell location

% % Extract dFF data from combinedTable, converting cell array to matrix
% dFF_data_sorted = cell2mat(combinedTable_sorted.dFF);
% % Apply z-score normalization to the dFF data
% dFF_data_zscore_sorted = zscore(dFF_data_sorted, [], 1);
% 
% %perivascular cells
% for perivascular_sorted = combinedTable_sorted.("Cell location (0,perivascular;1,adjacent;2,none)") == 0
%     dFF_data_Perivascular_sorted = combinedTable_sorted{perivascular_sorted,"dFF"};
% end
% dFF_data_Perivascular_sorted = cell2mat(dFF_data_Perivascular_sorted);
% 
% %Non-perivascular cells
% for non_perivascular_sorted = combinedTable_sorted.("Cell location (0,perivascular;1,adjacent;2,none)") == 2
%     dFF_data_NonPerivascular_sorted = combinedTable_sorted{non_perivascular_sorted,"dFF"};
% end
% dFF_data_NonPerivascular_sorted = cell2mat(dFF_data_NonPerivascular_sorted);
% 
% % Apply z-score normalization to the dFF data
% dFF_data_Perivascular_sorted_zscore = zscore(dFF_data_Perivascular_sorted, [], 1);
% dFF_data_NonPerivascular_sorted_zscore = zscore(dFF_data_NonPerivascular_sorted, [], 1);
% 
% % Combine the data from both groups
% combined_data_sorted = [dFF_data_Perivascular_sorted_zscore; dFF_data_NonPerivascular_sorted_zscore];
% 
% % Determine the colorbar limits based on the combined data
% caxis_limits_sorted = [min(combined_data_sorted(:)), max(combined_data_sorted(:))];
% 

%% Plots

% Heatmap all
figure;
imagesc(combined_data);
colormap(jet);
cb = colorbar; % Create colorbar
set(cb, 'Ticks', []); % Remove tick marks and numbers
caxis(caxis_limits_all_NM); % Set colorbar limits
xlabel('Time (sec)');
ylabel('Cell');

% Add a dark line under row 381
hold on;
line([0, size(combined_data, 2)], [381.5, 381.5], 'Color', 'k', 'LineWidth', 1);
hold off;

exportgraphics(gcf, 'heatmap3.eps', 'ContentType', 'vector');


figure;
imagesc(combined_data);
colormap(jet);
colorbar;
caxis(caxis_limits_all_NM); % Set colorbar limits
%title('dFF zscore');
xlabel('Time (sec)');
ylabel('Cell');
% Add a dark line under row 122
hold on;
line([0, size(combined_data, 2)], [381.5, 381.5], 'Color', 'k', 'LineWidth', 1);
hold off;

% Heatmap Perivascular
figure;
imagesc(dFF_data_Perivascular_zscore);
colormap(jet);
colorbar;
%caxis(caxis_limits_perivascular); % Set colorbar limits
caxis([minVal, maxVal]); % Set color axis
title('\DeltaF/F0 zscore Perivascular');
xlabel('Time(sec)');
ylabel('Cell');

% Heatmap NonPerivascular reduced number to match perivascular
figure;
imagesc(dFF_data_NonPerivascular_reduced_zscore);
colormap(jet);
colorbar;
%caxis(caxis_limits_nonPerivascular); % Set colorbar limits
caxis([minVal, maxVal]); % Set color axis
title('\DeltaF/F0 zscore NonPerivascular');
xlabel('Time(sec)');
ylabel('Cell');

% Heatmap NonPerivascular 
figure;
imagesc(dFF_data_NonPerivascular_zscore);
colormap(jet);
colorbar;
caxis(caxis_limits); % Set colorbar limits
title('dFF zscoreNonPerivascular');
xlabel('Time(sec)');
ylabel('Cell');

%% Order in number of events and plot dFF
% 
% % Sort combinedTable by number of events
% combinedTable_sorted = sortrows(combinedTable, 'Number of Events', 'descend');
% 
% %perivascular cells
% for perivascular_sorted = combinedTable_sorted.("Cell location (0,perivascular;1,adjacent;2,none)") == 0
%     dFF_data_Perivascular_sorted = combinedTable_sorted{perivascular_sorted,"dFF"};
% end
% dFF_data_Perivascular_sorted = cell2mat(dFF_data_Perivascular_sorted);
% 
% %Non-perivascular cells
% for non_perivascular_sorted = combinedTable_sorted.("Cell location (0,perivascular;1,adjacent;2,none)") == 2
%     dFF_data_NonPerivascular_sorted = combinedTable_sorted{non_perivascular_sorted,"dFF"};
% end
% dFF_data_NonPerivascular_sorted = cell2mat(dFF_data_NonPerivascular_sorted);
% dFF_data_NonPerivascular_sorted_reduced = dFF_data_NonPerivascular_sorted(1:122,:);
% 
% % Apply z-score normalization to the dFF data
% dFF_data_Perivascular_sorted_zscore = zscore(dFF_data_Perivascular_sorted, [], 1);
% dFF_data_NonPerivascular_sorted_zscore = zscore(dFF_data_NonPerivascular_sorted, [], 1);
% dFF_data_NonPerivascular_sorted_reduced_zscore = zscore(dFF_data_NonPerivascular_sorted_reduced, [], 1);
% 
% % Combine the data from both groups
% combined_data_sorted = [dFF_data_Perivascular_sorted_zscore; dFF_data_NonPerivascular_sorted_reduced_zscore];
% 
% % Determine the colorbar limits based on the combined data
% caxis_limits_sorted = [min(combined_data_sorted(:)), max(combined_data_sorted(:))];
% 
% % Heatmap Perivascular
% figure;
% imagesc(dFF_data_Perivascular_sorted_zscore);
% colormap(jet);
% colorbar;
% caxis(caxis_limits_sorted); % Set colorbar limits
% title('dF/F0 Perivascular');
% xlabel('Time(sec)');
% ylabel('Cell');
% 
% % Heatmap NonPerivascular
% figure;
% imagesc(dFF_data_NonPerivascular_sorted_reduced_zscore);
% colormap(jet);
% colorbar;
% caxis(caxis_limits_sorted); % Set colorbar limits
% title('dF/F0 NonPerivascular');
% xlabel('Time(sec)');
% ylabel('Cell');