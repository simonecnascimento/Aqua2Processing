function [adjMatrixCell, nodeDegreeINOUT, percentageConnected] = plotNodeConnectivityGraph(NetworkTable, data_CFU, data_analysis, simultaneousMatrixDelaybyCell_average, pwd, AquA_fileName)

    % Handle missing optional inputs
    if nargin < 6
        phase = []; % or default value (e.g., 1)
    end
    if nargin < 7
        tableNames = []; % or default value (e.g., {'Phase1'})
    end

    % Prepare adjacency matrix
    adjMatrixCell = cell2mat(simultaneousMatrixDelaybyCell_average);
    adjMatrixCell(isnan(adjMatrixCell)) = 0; % Replace NaNs with 0
    adjMatrixCell(1:size(adjMatrixCell,1)+1:end) = 0;

    % Remove double edges by keeping only i -> j when i < j
    n = size(adjMatrixCell, 1);
    for i = 1:n
        for j = i+1:n
            if adjMatrixCell(i, j) > 0 && adjMatrixCell(j, i) > 0
                % Remove j -> i
                adjMatrixCell(j, i) = 0;
            end
        end
    end

    % Create digraph
    G = digraph(adjMatrixCell);
    
    % Use only valid cells
    validCells = data_analysis.validCells;
    G.Nodes.Name = cellstr(string(validCells(:)));
    nodeNames = G.Nodes.Name;

    % Compute cell centers for valid cells
    centers_allCells = [];
    for idx = 1:numel(validCells)
        i = validCells(idx);
        cellMask = data_CFU.cfuInfo1{i,3};
        [rows1, cols1] = find(cellMask ~= 0);
        centers = [mean(rows1), mean(cols1)];
        centers_allCells = [centers_allCells; centers];
    end

    % Calculate out-degree for each node
    nodeDegree = outdegree(G);
    nodeDegreeINOUT = outdegree(G) + indegree(G);

    % Compute per-cell outgoing connectivity percentage within FOV
    nCells = NetworkTable.T.TotalCells;
    nodeDegreeINOUT = outdegree(G) + indegree(G);

    % remove CellsToExclude because these were from deleted events and multinucleated cells, which means they are invalid cells in the FOV
    % Get indices to remove from the table
    idxToRemove = NetworkTable.T.CellsToExclude{1}; 
    % Compute indices to keep
    idxToKeep = setdiff(1:length(nodeDegreeINOUT), idxToRemove);
    % Apply filtering
    nodeDegreeINOUT_clean = nodeDegreeINOUT(idxToKeep);
    percentageConnected = (nodeDegreeINOUT_clean / (nCells - 1)) * 100; % percentage connectivity per cell

    % Create figure
    digraphFig = figure;
    ax = axes; 
    set(ax, 'Color', 'w'); 
    hold on;

%     % in pixel
%     % Plot undirected graph with no arrows
%     h = plot(G, ...
%         'XData', centers_allCells(:,2), ...
%         'YData', centers_allCells(:,1), ...
%         'LineStyle', '-', ...
%         'Marker', 'o', ...
%         'NodeCData', nodeDegree, ...
%         'EdgeAlpha', 0.3, ...
%         'MarkerSize', 8, ...
%         'LineWidth', 1.5, ...
%         'NodeFontSize', 12, ...
%         'ArrowSize', 0.1); % Effectively hides arrows

    % in um
    pixelSize_um = 0.49;
    X_um = centers_allCells(:,2) * pixelSize_um;
    Y_um = centers_allCells(:,1) * pixelSize_um;

%     X_pixel = centers_allCells(:,2);
%     Y_pixel = centers_allCells(:,1);

     % Plot undirected graph with no arrows
    h = plot(G, ...
        'XData', X_um, ...
        'YData', Y_um, ...
        'LineStyle', '-', ...
        'Marker', 'o', ...
        'NodeCData', nodeDegreeINOUT, ...
        'EdgeAlpha', 0.3, ...
        'MarkerSize', 8, ...
        'LineWidth', 1.5, ...
        'NodeFontSize', 12, ...
        'ArrowSize', 0.1); % Effectively hides arrows

    colormap(ax, turbo);
    colorbar;
    title('Number of Connections');
    set(gca, 'YDir', 'reverse');
    xlim([1, 626 * pixelSize_um]);
    ylim([1, 422 * pixelSize_um]);

%     xlim([1, 626]);
%     ylim([1, 422]);
    hold off;

     % Save figure if savePath and AquA_fileName are provided
    if nargin > 3
        % Set up path
        fileTemp = extractBefore(AquA_fileName, "_AQuA2");
        pathTemp = extractBefore(pwd, "3.");
        
        % Define subfolder path
        subfolderDigraphName = 'figures\all cells (except multinucleated)\network_digraph\cell_digraphs\new';
        subfolderDigraphPath = fullfile(pathTemp, subfolderDigraphName);
        
        % Create the full file name with path
        if contains(pathTemp, 'CSD')
            digraphFilename = fullfile(subfolderDigraphPath, strcat(fileTemp, '_', tableNames{phase},'_diGraphCell.fig'));
        else
            digraphFilename = fullfile(subfolderDigraphPath, strcat(fileTemp,'_diGraphCell.fig'));
        end
        
        % Save the figure
        if ~exist(subfolderDigraphPath, 'dir')
            mkdir(subfolderDigraphPath); % Create the directory if it does not exist
        end
        saveas(digraphFig, digraphFilename);
    end

end




