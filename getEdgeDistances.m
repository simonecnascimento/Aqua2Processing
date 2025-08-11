function cellPairs_edges_distanceMicron_multipleAppearance = getEdgeDistances(adjMatrix, centers_allCells, numCells, rowsWithSingleAppearance)
    % Function to compute the pairwise distances between cells and add them
    % to the table with edge values from the adjacency matrix.
    
    % Initialize arrays for the values
    cellPairs_edges_distanceMicron = [];
    
    % Initialize distance matrix in pixels (for later conversion)
    distancesPixels = zeros(numCells, numCells);
    
    % Loop through all pairs of nodes
    for i = 1:numCells
        % Calculate the pairwise distances
        for j = i:numCells
            % Calculate the Euclidean distance between nodes i and j (in pixels)
            distancesPixels(i,j) = sqrt((centers_allCells(i,1) - centers_allCells(j,1))^2 + (centers_allCells(i,2) - centers_allCells(j,2))^2);
            distancesPixels(j,i) = distancesPixels(i,j);  % The matrix is symmetric
        end
    end
    
    % Convert pixels to microns
    distancesMicrons = distancesPixels * 0.49;

    % Loop through each cell pair and add distances where adjacency is non-zero
    for i = 1:numCells
        for j = 1:numCells
            if adjMatrix(i, j) ~= 0 % Ensure it's not zero (valid edge)
                % Prepare the edge data
                cellPair = [i, j];
                edge_value = adjMatrix(i, j);  % Get the edge value from the adjacency matrix
                
                % Calculate the distance for the pair
                distanceCellPair = distancesMicrons(i, j);
                
                % Create the edge data including distance
                cellPair_edge_distanceMicron = [cellPair, edge_value, distanceCellPair];
                
                % Store the result in the final output array
                cellPairs_edges_distanceMicron = [cellPairs_edges_distanceMicron; cellPair_edge_distanceMicron];
            end
        end
    end

    cellPairs_edges_distanceMicron_multipleAppearance = cellPairs_edges_distanceMicron;
    %cellPairs_edges_distanceMicron_multipleAppearance(rowsWithSingleAppearance, :) = [];
end


% 
% 
% % Initialize arrays for the values
% edge_values = [];
% cellPairs_edges = [];
% cellPair_edge_distanceMicron = [];
% cellPairs_edges_distanceMicron = [];
% 
% % Loop through all pairs of nodes
% for i = 1:numCells
%     % Calculate the pairwise distances
%     for j = i:numCells
%         % Calculate the Euclidean distance between nodes i and j
%         distancesPixels(i,j) = sqrt((centers_allCells(i,1) - centers_allCells(j,1))^2 + (centers_allCells(i,2) - centers_allCells(j,2))^2);
%         distancesPixels(j,i) = distancesPixels(i,j);  % The matrix is symmetric
%         % convert pixel to um
%         distancesMicrons = distancesPixels * 0.49;
%     end
% 
%     for j = 1:numCells 
%         if adjMatrix(i, j) ~= 0 % Ensure it's not zero
%             cellPair = [i,j];
%             edge_values = adjMatrix(i, j);
%             cellPair_edge = [cellPair, edge_values];
%             cellPairs_edges = [cellPairs_edges; cellPair_edge];
% 
%             distanceCellPair = distancesMicrons(i, j);
%             cellPair_edge_distanceMicron = [cellPair_edge, distanceCellPair];
%             cellPairs_edges_distanceMicron = [cellPairs_edges_distanceMicron; cellPair_edge_distanceMicron];
% 
%         end
%     end
% 
% end


