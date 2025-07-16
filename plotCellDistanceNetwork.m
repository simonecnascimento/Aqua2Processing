function [adjMatrix, centers_allCells, G, rowsWithSingleAppearance, nodeDegree] = plotCellDistanceNetwork(numCells, data_CFU, data_analysis, simultaneousMatrixDelaybyCell, simultaneousMatrixDelaybyCell_average, pwd, AquA_fileName, tableNames, phase)
    % Function to plot a directed graph representing cell distance networks
    % and optionally save the resulting figure.
    % Inputs:
    % - data_CFU: Data containing cell coordinates
    % - data_analysis: Analysis data including perivascularCells
    % - simultaneousMatrixDelaybyCell2: Cell array with edge weights
    % - simultaneousMatrixDelaybyCell_average2: Cell array with adjacency matrix values


    % ---- Preprocessing ----
    simultaneousMatrixDelaybyCell_average2 = simultaneousMatrixDelaybyCell_average;  
    simultaneousMatrixDelaybyCell2 = simultaneousMatrixDelaybyCell;
    
    % Process each cell to remove zeros
    for cellX = 1:size(simultaneousMatrixDelaybyCell2) %rows
        for cellY = 1:size(simultaneousMatrixDelaybyCell2) %columns
            currentCell = simultaneousMatrixDelaybyCell2{cellX, cellY};
            if isnumeric(currentCell) % Ensure the cell contains numeric data
                % Remove zeros and keep non-zero elements
                simultaneousMatrixDelaybyCell2{cellX, cellY} = currentCell(currentCell ~= 0);
            end
            if cellX == cellY
                simultaneousMatrixDelaybyCell2{cellX, cellY} = [];
                simultaneousMatrixDelaybyCell_average2{cellX, cellY} = NaN;
            end
        end
    end
    % ---- End Preprocessing ----

    % Extract cell coordinates
    cellCoords = [data_CFU.cfuInfo1(:,1), data_CFU.cfuInfo1(:,3)];
    numCells = size(cellCoords, 1);
    
    % Compute cell centers
    centers_allCells = [];
    for i = 1:numCells
        [rows1, cols1] = find(cellCoords{i,2} ~= 0); % Extract coordinates of non-zero elements
        centers = [mean(rows1), mean(cols1)]; % Compute center coordinates
        centers_allCells = [centers_allCells; centers];
    end
    
    % Nodes, edges, and weights
    edgeWeights = simultaneousMatrixDelaybyCell2;
    connectionsMatrix = simultaneousMatrixDelaybyCell_average2;

    % Calculate edge weights based on the number of elements in each cell
    numericWeights = cellfun(@(c) numel(c), edgeWeights);

    % Process adjacency matrix
    adjMatrix = cell2mat(connectionsMatrix);
    adjMatrix(isnan(adjMatrix)) = 0; % Replace NaN with 0
    %adjMatrix(numericWeights < 2) = 0; % Remove edges with weights < 2

    % Create directed graph
    G = digraph(adjMatrix);
    
    % Get node degree
    nodeDegree = outdegree(G);

    % Identify nodes with single appearance
    endNodes = G.Edges.EndNodes;
    allNodes = endNodes(:);

    % Check if allNodes is empty and skip node processing if true
    if isempty(allNodes)
        % Set up figure without processing nodes
        disp('allNodes is empty, skipping node processing.');
    else
        nodeCounts = histcounts(allNodes, 1:(max(allNodes)+1));
        nodesWithSingleAppearance = find(nodeCounts == 1);
    
        % Find indices of edges connected to nodes with single appearance
        rowsWithSingleAppearance = [];
        for i = 1:length(nodesWithSingleAppearance)
            node = nodesWithSingleAppearance(i);
            indices = find(endNodes(:, 1) == node | endNodes(:, 2) == node);
            rowsWithSingleAppearance = [rowsWithSingleAppearance; indices];
        end
    end

    % Set up figure
    digraphFig = figure;
    axesHandle = axes;
    set(axesHandle, 'Color', 'w');
    hold(axesHandle, 'on');
    
    % Plot graph
    h = plot(G, 'XData', centers_allCells(:, 2), 'YData', centers_allCells(:, 1));
    h.MarkerSize = 7;
    h.ArrowSize = 10;
    h.NodeFontSize = 15;
    h.EdgeColor = [0.4, 0.4, 0.4];

    % Assign node colors based on degree
    h.NodeCData = nodeDegree; % NodeCData will color nodes by their degree
    colormap(axesHandle, spring); % Colormap for degree visualization
    colorbar; % Add a colorbar for reference

%     % Set line thickness for edges
%     validWeights = numericWeights(numericWeights >= 2);
%     h.LineWidth = validWeights * 2;

    % Assign node shape based on type of cell (perivascular vs non-perivascular)
    perivascular = data_analysis.perivascularCells;
    nonPerivascular = setdiff(1:numCells, perivascular); % Other nodes

    % Highlight perivascular nodes (stars)
    highlight(h, perivascular, 'Marker', 'p', 'MarkerSize', 10);  
    % Highlight non-perivascular nodes (round markers)
    highlight(h, nonPerivascular, 'Marker', 'o', 'MarkerSize', 7); 

%     nodeTypes = ismember(1:numCells, perivascular);
%     nodeColors = 2 * (~nodeTypes); % Red (0) or Blue (2)
%     colormap(axesHandle, [1 0 0; 0 0 1]);
%     h.NodeCData = nodeColors;

    % Adjust axis
    title('Cell Distance Network on Image Frame');
    set(gca, 'YDir', 'reverse'); % 'reverse'
    xlim([1, 626]);
    ylim([1, 422]);

    % Annotate edges with weights and style
    for i = 1:numedges(G)
        source = G.Edges.EndNodes(i, 1);
        target = G.Edges.EndNodes(i, 2);
        sourceCoord = h.XData(source);
        targetCoord = h.XData(target);
        ySourceCoord = h.YData(source);
        yTargetCoord = h.YData(target);
        midpointX = (sourceCoord + targetCoord) / 2;
        midpointY = (ySourceCoord + yTargetCoord) / 2;
        edgeWeight = G.Edges.Weight(i);
        edgeWeightInt = round(edgeWeight);

        % Check if edge is dashed
        if ismember(i, rowsWithSingleAppearance)
            line([sourceCoord, targetCoord], [ySourceCoord, yTargetCoord], ...
                'LineStyle', '--', 'LineWidth', 2, 'Color', [0.4, 0.4, 0.4]);
        end

        % Display weight
%         text(midpointX, midpointY, num2str(edgeWeightInt), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
%             'FontSize', 12, 'Color', 'k');
    end
    hold off;

    % Save figure if savePath and AquA_fileName are provided
    if nargin > 4
        % Set up path
        fileTemp = extractBefore(AquA_fileName, "_AQuA2");
        pathTemp = extractBefore(pwd, "3.");
        
        % Define subfolder path
        subfolderDigraphName = 'figures\all cells (except multinucleated)\network_digraph\digraphs\new';
        subfolderDigraphPath = fullfile(pathTemp, subfolderDigraphName);
        
        % Create the full file name with path
        if contains(pathTemp, 'CSD')
            digraphFilename = fullfile(subfolderDigraphPath, strcat(fileTemp, '_', tableNames{phase},'_diGraph.fig'));
        else
            digraphFilename = fullfile(subfolderDigraphPath, strcat(fileTemp,'_diGraph.fig'));
        end
        
        % Save the figure
        if ~exist(subfolderDigraphPath, 'dir')
            mkdir(subfolderDigraphPath); % Create the directory if it does not exist
        end
        saveas(digraphFig, digraphFilename);
    end
end
