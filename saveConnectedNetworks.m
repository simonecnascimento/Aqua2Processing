function saveConnectedNetworks(G, centers_allCells, pwd, AquA_fileName, tableNames, phase)
    % SAVECONNECTEDNETWORKS Identifies and visualizes disconnected networks in a graph.
    %
    % Inputs:
    %   G               - A directed graph object (digraph).
    %   pwd             - Current directory path.
    %   AquA_fileName   - File name for saving the figure.

    % Create a figure for connected networks
    connected_network = figure;
    h = plot(G, 'XData', centers_allCells(:, 2), 'YData', centers_allCells(:, 1));
    set(gca, 'YDir', 'reverse');

    % Identify connected components - check for multiple disconnected networks
    [components, numComponents] = conncomp(G, 'Type', 'strong');

    % Ensure numComponents is a scalar
    numComponents = max(components);

    % Assign a unique color to each component
    h.NodeCData = components;

    % Set colormap and add colorbar
    colormap(jet(numComponents));
    colorbar;

    % Set title
    title('Strongly Connected Components');

    % Generate paths for saving the figure
    fileTemp = extractBefore(AquA_fileName, "_AQuA2");
    pathTemp = extractBefore(pwd, "3.");

    % Define subfolder path
    subfolderNetworkName = 'figures\all cells (except multinucleated)\network_digraph\networks\new';
    subfolderNetworkPath = fullfile(pathTemp, subfolderNetworkName);

    % Create the full file name with path
    if contains(pathTemp, 'CSD')
       networksFilename = fullfile(subfolderNetworkPath, strcat(fileTemp, '_', tableNames{phase}, '_networks.fig'));
    else
       networksFilename = fullfile(subfolderNetworkPath, strcat(fileTemp, '_networks.fig'));
    end

    % Ensure the subfolder exists, create it if necessary
    if ~exist(subfolderNetworkPath, 'dir')
        mkdir(subfolderNetworkPath);
    end

    % Save the figure
    saveas(connected_network, networksFilename);
end
