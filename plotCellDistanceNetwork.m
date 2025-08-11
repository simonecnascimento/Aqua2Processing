function [adjMatrix, centers_allCells, G, rowsWithSingleAppearance, nodeDegree] = plotCellDistanceNetwork(data_CFU, data_analysis, simultaneousMatrixDelaybyCell, simultaneousMatrixDelaybyCell_average, pwd, AquA_fileName, tableNames, phase)
    % Function to plot a directed graph representing cell distance networks
    % and optionally save the resulting figure.
    % Inputs:
    % - data_CFU: Data containing cell coordinates
    % - data_analysis: Analysis data including perivascularCells
    % - simultaneousMatrixDelaybyCell2: Cell array with edge weights
    % - simultaneousMatrixDelaybyCell_average2: Cell array with adjacency matrix values


    % Handle missing optional inputs
    if nargin < 8
        phase = []; % or default value (e.g., 1)
    end
    if nargin < 7
        tableNames = []; % or default value (e.g., {'Phase1'})
    end

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

    % Select only the valid cells
    validCells = data_analysis.validCells;
    
    % Compute cell centers
    centers_allCells = [];
    for idx = 1:numel(validCells)
        i = validCells(idx);  
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
    
    % Use only valid cells
    G.Nodes.Name = cellstr(string(validCells(:)));
    nodeNames = G.Nodes.Name;

    % Get node degree
    nodeDegree = outdegree(G);

    % Identify nodes with single appearance
    endNodes = G.Edges.EndNodes;

    % Convert to numeric arrays
    sourceNodes = str2double(endNodes(:,1));
    targetNodes = str2double(endNodes(:,2)); 
    % Combine all node indices
    allNodes = [sourceNodes; targetNodes];

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
            indices = find(sourceNodes == node | targetNodes == node);
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

    % Convert numeric node indices to strings
    nonPerivascularNames = arrayfun(@num2str, nonPerivascular, 'UniformOutput', false);
    perivascularNames = arrayfun(@num2str, perivascular, 'UniformOutput', false);

    % Find valid nodes in graph for each group
    validNonPerivascular = intersect(nonPerivascularNames, nodeNames);
    validPerivascular = intersect(perivascularNames, nodeNames);
    
    % Highlight nonPerivascular nodes
    highlight(h, validNonPerivascular, 'Marker', 'o', 'MarkerSize', 7);
    
    % Highlight perivascular nodes
    highlight(h, validPerivascular, 'Marker', 's', 'MarkerSize', 8);

%     % Highlight perivascular nodes (stars)
%     highlight(h, nodeNames(perivascular), 'Marker', 'p', 'MarkerSize', 10);  
%     % Highlight non-perivascular nodes (round markers)
%     highlight(h, nodeNames(nonPerivascular), 'Marker', 'o', 'MarkerSize', 7);
%     %highlight(h, nonPerivascular, 'Marker', 'o', 'MarkerSize', 7); 

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
        sourceName = G.Edges.EndNodes{i, 1}; % e.g. '10'
        targetName = G.Edges.EndNodes{i, 2}; % e.g. '5'
        
        % Find index of sourceName in G.Nodes.Name
        [~, sourceIdx] = ismember(sourceName, G.Nodes.Name);
        [~, targetIdx] = ismember(targetName, G.Nodes.Name);
        
        % Now use numeric indices to get coordinates
        sourceCoordX = h.XData(sourceIdx);
        sourceCoordY = h.YData(sourceIdx);
        
        targetCoordX = h.XData(targetIdx);
        targetCoordY = h.YData(targetIdx);

        midpointX = (sourceCoordX + targetCoordX) / 2;
        midpointY = (sourceCoordY + targetCoordY) / 2;
        edgeWeight = G.Edges.Weight(i);
        edgeWeightInt = round(edgeWeight);

        % Check if edge is dashed
        if ismember(i, rowsWithSingleAppearance)
            line([sourceCoordX, targetCoordX], [sourceCoordY, targetCoordY], ...
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
        subfolderDigraphName = 'figures\all cells (except multinucleated)\network_digraph\event_digraphs\new';
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
