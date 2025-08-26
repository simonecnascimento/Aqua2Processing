function rasterTable_sorted = plotRasterByCluster(rasterTable, selectedClusters, experiment)
% plotRasterByCluster - Creates a raster plot from preCSD, duringCSD, and postCSD event times.
% Only plots cells belonging to specified clusters.
%
% Syntax:
%   rasterTable_sorted = plotRasterByCluster(rasterTable, selectedClusters)
%
% Inputs:
%   rasterTable      - A table with columns:
%                      [cellID, preCSD, duringCSD, postCSD, baseline_preCSD, clusterID]
%   selectedClusters - Vector of cluster IDs to include (e.g. [2,4,9])

    % Extract and sort clusters
    
    %Cluster 2 (acute increase)
    cluster2_only = rasterTable(rasterTable.eventRate_clusterID == 2, :);
    % Get number of elements in each time period
    nDuring = cellfun(@numel, cluster2_only.duringCSD);
    nPre    = cellfun(@numel, cluster2_only.preCSD);
    nPost   = cellfun(@numel, cluster2_only.postCSD);
    % Sort by duringCSD, then preCSD, then postCSD (all descending)
    [~, sortIdx] = sortrows([nDuring, nPre, nPost], [-1, -2, -3]);
    cluster2_sorted = cluster2_only(sortIdx, :);

%     if strcmp(experiment, 'CSD')
%         cluster2_sorted(27:28, :) = []; %CSD
%     elseif strcmp(experiment, 'BIBN-CSD')
%         cluster2_sorted(1:3,:) = []; %BIBN
%     end

    %Cluster 4 (chronic increase)
    cluster4_only = rasterTable(rasterTable.eventRate_clusterID == 4, :);
    % Compute number of elements in postCSD and preCSD
    nPost = cellfun(@numel, cluster4_only.postCSD);
    nPre  = cellfun(@numel, cluster4_only.preCSD);
    % Sort first by postCSD, then by preCSD (both descending)
    [~, sortIdx] = sortrows([nPost, nPre], [-1, -2]);
    % Apply sorting
    cluster4_sorted = cluster4_only(sortIdx, :);
    
    % Cluster6
    cluster6_only = rasterTable(rasterTable.eventRate_clusterID == 6, :);
    lens_pre = cellfun(@numel, cluster6_only.preCSD);
    lens_post = cellfun(@numel, cluster6_only.postCSD);
    [~, sortIdx] = sortrows([lens_pre(:), lens_post(:)], [-1, -2]); % Descending
    cluster6_sorted = cluster6_only(sortIdx, :);
    % Shuffle cluster 6 subset
    if strcmp(experiment, 'CSD')
        rangeToShuffle = 87:135; % CSD
    elseif strcmp(experiment, 'BIBN-CSD')
        rangeToShuffle = 14:25; % BIBN
    end
    shuffledIndices = rangeToShuffle(randperm(length(rangeToShuffle)));
    cluster6_sorted(rangeToShuffle, :) = cluster6_sorted(shuffledIndices, :);

%     % Cluster9
%     cluster9_only = rasterTable(rasterTable.eventRate_clusterID == 9, :);
%     [~, sortIdx] = sort(cellfun(@numel, cluster9_only.preCSD), 'descend');
%     cluster9_sorted = cluster9_only(sortIdx, :);
% 
%     % Shuffle cluster 9 subset
%     rangeToShuffle = 76:100; % CSD
%     rangeToShuffle = 14:25; % BIBN
%     shuffledIndices = rangeToShuffle(randperm(length(rangeToShuffle)));
%     cluster9_sorted(rangeToShuffle, :) = cluster9_sorted(shuffledIndices, :);

    % Concatenate
    rasterTable_sorted = [cluster2_sorted; cluster4_sorted; cluster6_sorted]; %cluster4_sorted; 

    if strcmp(experiment, 'CSD')
        rasterTable_sorted([1:2,4:7],:) = []; %CSD
    elseif strcmp(experiment, 'BIBN-CSD')
        rasterTable_sorted(1:2,:) = []; %BIBN
    end



    % Initialize figure
    figure;
    hold on;

    % Determine number of rows
    rowCounter = height(rasterTable_sorted);

    % --- Draw gray patch BEFORE lines (so it goes behind)
    grayX = [1854.5, 1917.5, 1917.5, 1854.5];
    grayY = [0, 0, rowCounter, rowCounter];
    patch(grayX, grayY, [0.85 0.85 0.85], 'EdgeColor', 'none');

    % --- Plot raster lines
    for i = 1:rowCounter
        times = [rasterTable_sorted.preCSD{i}, ...
                 rasterTable_sorted.duringCSD{i}, ...
                 rasterTable_sorted.postCSD{i}];

        for t = times
            if t <= 3772
                line([t t], [i - 0.4, i + 0.4], 'Color', 'k', 'LineWidth', 1);
            end
        end
    end

    % --- Cluster divider lines
    clusterOrder = rasterTable_sorted.eventRate_clusterID;
    clusterBoundaries = find(diff(clusterOrder));
    for i = 1:length(clusterBoundaries)
        yline(clusterBoundaries(i)+0.5, 'Color', 'r', 'LineWidth', 0.5);
    end

    % --- Format axes (Time in minutes)
    xlabel('Time (min)', 'FontSize', 6);
    ylabel('Cell index', 'FontSize', 6);
    set(gca, 'FontSize', 6);

    ylim([0 rowCounter+1]);
    xlim([-20 3790]);

    % --- Set ticks in minutes (10 min intervals)
    framePerMin = 927 / 15; % ≈ 61.8
    tickFrames = 0:(10 * framePerMin):3772;
    tickMins = round(tickFrames / framePerMin);
    set(gca, 'XTick', tickFrames, 'XTickLabel', tickMins);

    % --- Axis direction and appearance
    set(gca, 'YDir', 'reverse', 'TickLength', [0 0]);
    box off;

    % --- Title
    title(['Raster Plot (pre/during/post CSD) for Clusters: ', num2str(selectedClusters)], 'FontSize', 6);

    % --- Optional: Export for Illustrator
    print(gcf, 'raster_plot_BIBN', '-depsc', '-painters');

end