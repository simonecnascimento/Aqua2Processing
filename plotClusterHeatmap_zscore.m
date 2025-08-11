function plotClusterHeatmap_zscore(combinedTable, includedClusters, savePath)
% plotClusterHeatmap - Plots a heatmap of calcium traces grouped by cluster
%
% Syntax:
%   plotClusterHeatmap(combinedTable, clusterCol, traceCol, includedClusters)
%
% Inputs:
%   combinedTable      - MATLAB table containing all cell data
%   includedClusters   - Vector of cluster IDs to include (e.g. [2,4,8,9])

   
    % Initialize
    selectedData = {};
    selectedClusters = [];

    clusterCol = combinedTable{:,17};
    traceCol = combinedTable{:,"dFF"};
    filenameCol = combinedTable{:,"fileNameColumn"};

    % Loop through each cluster
    for i = 1:length(includedClusters)
        clust = includedClusters(i);
        idx = find(clusterCol == clust);

%         % Limit to max 10 cells per cluster
%         if numel(idx) >  10
%             idx = idx(1:10); %  idx = idx(1:10); or use randperm(numel(idx), 10) for random
%         end

        % Append
        selectedData = [selectedData; traceCol(idx)];
        selectedClusters = [selectedClusters; repmat(clust, length(idx), 1)];
    end

    % Convert to matrix
    dFF_Matrix = cell2mat(selectedData);

    for x = 1:size(dFF_Matrix,1)
        figure;
        plot(dFF_Matrix(x,:))
    end

    toRemove_4 = any(dFF_Matrix > 4, 2);  % Find rows with any value > 4
    dFF_Matrix(toRemove_4, :) = [];       
    selectedClusters(toRemove_4, :) = [];  
    selectedData(toRemove_4, :) = [];      

    toRemove_3 = any(dFF_Matrix > 3, 2);  % Find rows with any value > 3
    dFF_Matrix(toRemove_3, :) = [];     

    toRemove = 64:80;
    dFF_Matrix(toRemove, :) = [];
    selectedClusters(toRemove, :) = [];  
    selectedData(toRemove, :) = [];   


%     % Z-score by row (each cell individually)
%     zscored_dFF = zscore(dFF_Matrix, 0, 2);  % 0 = normalize using std, 2 = operate along rows

%     % Z-score by median "This signal was standardized akin to a z-score
%     % operation by subtracting the median value and dividing by the standard deviation 
%     % (calculated during quiet wakefulness, an epoch with low levels of evoked activity)."
%     
%     baselineFrames = 1:1854;  % Change this to match your quiet epoch
%     % Preallocate
%     zscored_dFF = zeros(size(dFF_Matrix));
%     
%     % Loop over each cell
%     for i = 1:size(dFF_Matrix, 1)
%         trace = dFF_Matrix(i, :);
%         baselineStd = std(trace(baselineFrames));  % std from quiet period
%         baselineMedian = median(trace);         % median from full trace
%     
%         zscored_dFF(i, :) = (trace - baselineMedian) / baselineStd;
%     end

    %remove = [6,7,46];
    %zscored_dFF(remove,:) = []; 

%     [maxVals, maxIndices] = max(zscored_dFF, [], 2);  % For each row (cell), get max and its index
%     [minVals, minIndices] = min(zscored_dFF, [], 2);  % For each row (cell), get max and its index
%     histMax = figure;
%     histogram(maxIndices,60)
%     histMin = figure;
%     histogram(minIndices,60)

    % Plot
    figure;
    imagesc(dFF_Matrix(:,1:3772))
    colormap(flipud(gray));
    colorbar;
    xlabel('Time (frames)');
    %ylabel('Cells (sorted by cluster)');
    title(['Heatmap clusters: ', num2str(includedClusters)]);

    % Cluster dividers
    hold on;
    clusterBoundaries = find(diff(selectedClusters));
    for i = 1:length(clusterBoundaries)
        yline(clusterBoundaries(i)+0.5, 'Color', 'k', 'LineWidth', 2);
    end
    
    % Define specific transitions you want to mark
    targetTransitions = [3 1; 7 3]; % CSD [3 1; 7 3] [2 1; 7 3] [3 1; 7 3]
    %targetTransitions = [2 1; 7 6]; % BIBN-CSD
    
    % Loop through and add lines where these transitions occur
    for i = 1:size(selectedClusters, 1)-1
        current = selectedClusters(i);
        next = selectedClusters(i+1);
        
        for j = 1:size(targetTransitions, 1)
            if current == targetTransitions(j, 1) && next == targetTransitions(j, 2)
                yline(i+0.5, 'w', 'LineWidth', 2);  % red line for visibility
            end
        end
    end


    % Vertical lines at frames 1854 and 1855
    xline(1854, '--r', 'LineWidth', 1.5);

    % Optional save
    if nargin >= 5 && ~isempty(savePath)
        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end
        filename = fullfile(savePath, 'clusterHeatmap.eps');
        saveas(gcf, filename);
        fprintf('Figure saved to: %s\n', filename);
    end
end

