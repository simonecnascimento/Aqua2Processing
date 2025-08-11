function resultantVector = sumEdgeVectors(G, nodePositions, rowsWithSingleAppearance, pwd, AquA_fileName, tableNames, phase)
    % SUMEDGEVECTORS Sums all edge vectors in a directed graph and computes the resultant vector.
    %
    % Input:
    %   G - A directed graph object (digraph)
    %   nodePositions - An N x 2 matrix of [x, y] positions for each node
    %                   where N is the number of nodes in G.
    %
    % Output:
    %   resultantVector - The resultant vector [x, y] after summing all edge vectors.

    % Ensure input is a digraph
    if ~isa(G, 'digraph')
        error('Input must be a digraph object.');
    end

    % Ensure nodePositions match the number of nodes
    if size(nodePositions, 1) ~= numnodes(G)
        error('nodePositions must have one row per node in G.');
    end

    % Initialize the resultant vector
    resultantVector = [0, 0];

    % Loop through each edge and compute its vector
    for i = 1:height(G.Edges)
        if ismember(i, rowsWithSingleAppearance)
            resultantVector = resultantVector;
        else
            % Get source and target node indices
            sourceNode = findnode(G, G.Edges.EndNodes(i, 1));
            targetNode = findnode(G, G.Edges.EndNodes(i, 2));
    
            % Get source and target positions
            sourcePos = nodePositions(sourceNode, :);
            targetPos = nodePositions(targetNode, :);
    
            % Compute the vector for this edge
            edgeVector = sourcePos - targetPos;
            modifiedEdgeVector = [edgeVector(1), -edgeVector(2)];
    
            % Add it to the resultant vector
            resultantVector = resultantVector + modifiedEdgeVector;
        end
    end

    % Display the resultant vector
    disp('Resultant Vector:');
    disp(resultantVector);

    % Plot the resultant vector
    vectorFig = figure;
    hold on;
    % Plot origin
    plot(0, 0, 'ro', 'MarkerSize', 10, 'DisplayName', 'Origin');
    % Plot the resultant vector as an arrow
    quiver(0, 0, resultantVector(2), resultantVector(1), ...
        0, 'MaxHeadSize', 0.5, 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'Resultant Vector');
    axis equal;
    xlabel('X');
    ylabel('Y');
    title('Resultant Vector from Summing Edge Vectors');
    legend;
    hold off;

    % Save the resultant vector plot if a save path is provided
    if nargin > 3
        % Set up path
        fileTemp = extractBefore(AquA_fileName, "_AQuA2");
        pathTemp = extractBefore(pwd, "3.");
        
        % Define subfolder path
        subfolderVectorName = 'figures\all cells (except multinucleated)\network_digraph\resultingVector\new';
        subfolderVectorPath = fullfile(pathTemp, subfolderVectorName);
        
        % Create the full file name with path
        if contains(pathTemp, 'CSD')
            vectorFilename = fullfile(subfolderVectorPath, strcat(fileTemp, '_', tableNames{phase}, '_resultingVector.fig'));
        else
            vectorFilename = fullfile(subfolderVectorPath, strcat(fileTemp, '_resultingVector.fig'));
        end
        
        % Save the figure
        if ~exist(subfolderVectorPath, 'dir')
            mkdir(subfolderVectorPath); % Create the directory if it does not exist
        end
        saveas(vectorFig, vectorFilename);
    end

end
