%% Load .mat files 

clear all;

% Set the directory for the experiment you need
fullCraniotomyBaselineDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\3._analysisByEvent.mat';
%fullCraniotomyCSDDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick\1._analysis.mat';
% thinBoneBaselineDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\thinBone\_analysisByEvent.mat';

cd (fullCraniotomyBaselineDir);
% Get all .mat files in the directory
FilesAll = dir(fullfile(fullCraniotomyBaselineDir, '*_analysisByEvent.mat')); 

% Extract file names
fileNames = {FilesAll.name};

% Initialize an array to store the extracted parts for sorting
sortKey = [];

% Loop through all filenames to extract the parts for sorting
for file = 1:length(fileNames)
    filename = fileNames{file};
    
    % Extract the number after "Pf4Ai162-" (e.g., '2' from 'Pf4Ai162-2')
    numberAfterPrefix = sscanf(filename, 'Pf4Ai162-%d', 1);
    
    % Extract the date (e.g., '221130' from 'Pf4Ai162-2_221130_FOV6')
    dateStr = regexp(filename, '\d{6}', 'match', 'once');
    
    % Extract the FOV number (e.g., 'FOV6' from 'Pf4Ai162-2_221130_FOV6')
    fovNumber = sscanf(filename, 'Pf4Ai162-%*d_%*d_FOV%d', 1);
    
    % Store the extracted values in a matrix for sorting
    sortKey = [sortKey; numberAfterPrefix, str2double(dateStr), fovNumber];
end

% Sort by the three columns: numberAfterPrefix, date, fovNumber
[~, idx] = sortrows(sortKey);
sortedFileNames = fileNames(idx);

% Initialize an empty table
networkTable_all = table();

%% extract network spatial density

% Initialize a cell array to store the upperTriAverages for each experiment, for distance between cells

allUpperTriValues = cell(length(sortedFileNames), 1); %distance between pair of cells
eventDuration_simultaneousEvents_all = cell(length(sortedFileNames), 1);
simultaneousMatrixDelaybyCell_average_all = cell(length(sortedFileNames), 1);
listOfEvents_perCell_nodeOUTdegree_all = cell(length(sortedFileNames), 1);
cellPairs_edges_distanceMicron_multipleAppearance_all = [];

for experiment = 1:length(sortedFileNames)
 
    [data_analysis, data_aqua, data_CFU, AquA_fileName] = loadAnalysisData(sortedFileNames, experiment);

    % FILL TABLE WITH CELL NUMBERS, MAPS, AND EVENTS PROPAGATIONS
    % get network data from aqua
    nSimultanousEvents = num2cell(data_aqua.res.fts1.networkAll.nOccurSameTime);
    simultaneousEvents_all = data_aqua.res.fts1.networkAll.occurSameTimeList;
    numEvents = height(simultaneousEvents_all); % Define the number of events

    % for current event
    propagationMap_event = cell(numEvents,1);
    % for cell corresponding the current event
    cellNumber_event = cell(numEvents,1);
    cellType_event = cell(numEvents,1);
    cellMap_event = cell(numEvents,1);
    
    % for all simultaneous events
    propagationMap_all = cell(numEvents,1);
    % for cells corresponding the all simultaneous events
    cellNumber_all = cell(numEvents,1);
    cellType_all = cell(numEvents,1);
    cellMap_all = cell(numEvents,1);
    % for cells corresponding the all simultaneous events, but without repeating cell numbers
    cellNumber_all_unique = cell(numEvents,1);
    cellMap_all_unique = cell(numEvents,1);

    % for events removing current Event
    simultaneousEvents_network = simultaneousEvents_all; %duplicate 
    simultaneousEvents_network_corrected = simultaneousEvents_all; %duplicate and remove any eventToDelete
    propagationMap_network = cell(numEvents,1);
    % for cell corresponding to events in the network
    cellNumber_network = cell(numEvents,1);
    cellType_network = cell(numEvents,1);
    cellMap_network = cell(numEvents,1);

    % combine propagation map of event to map of cell
    eventNetwork_cellMap_all = cell(numEvents,1);
    % calculate shortest distance between centers
    shortestDistanceBetweenCenters = cell(numEvents,1);

    % create a table with network results
    networkData = [propagationMap_event, cellNumber_event, cellType_event, cellMap_event, nSimultanousEvents, ...
        simultaneousEvents_all, propagationMap_all, cellNumber_all, cellType_all, cellMap_all, ...
        cellNumber_all_unique, cellMap_all_unique, eventNetwork_cellMap_all, shortestDistanceBetweenCenters, ...
        simultaneousEvents_network, simultaneousEvents_network_corrected, propagationMap_network, cellNumber_network, cellType_network, cellMap_network];
    networkData = cell2table(networkData);
    networkData.Properties.VariableNames = {'propagationMap_event','cellNumber_event', 'cellType_event','cellMap_event', 'nSimultanousEvents', ...
        'simultaneousEvents_all', 'propagationMap_all', 'cellNumber_all', 'cellType_all', 'cellMap_all', 'cellNumber_all_unique', 'cellMap_all_unique', 'eventNetwork_cellMap_all', 'shortestDistanceBetweenCenters', ...
        'simultaneousEvents_network', 'simultaneousEvents_network_corrected', 'propagationMap_network', 'cellNumber_network', 'cellType_network', 'cellMap_network'};
      
%     % Duplicate table to remove repetitive cells in the future
%     networkData_cleaned = networkData;
% %     cellNumbers = num2cell(networkData_cleaned{currentEvent, 'cellNumber_all'}{1}); % Assuming cellNumber_all is a cell array
% %     cellMaps = networkData_cleaned{currentEvent, 'cellMap_all'}{1}; % cellMap_all might already be a cell array

    % Matrix by event
    simultaneousMatrixByEvent = zeros(numEvents, numEvents);
    [numRows, numCols] = size(simultaneousMatrixByEvent);
    simultaneousMatrixDelayByEvent = zeros(numEvents, numEvents);
    startingFrame = data_analysis.resultsRaw.ftsTb(3,:); %starting frame of event to compute delay
    
    % Matrix by cell
    numCells = size(data_CFU.cfuInfo1,1);
    simultaneousMatrixDelaybyCell = cell(numCells, numCells);
    simultaneousMatrixDelaybyCell_average = cell(numCells, numCells);

    for currentEvent = 1:size(networkData,1)

        % find the propagation matrix related to the current event
        propagationMap_event = data_aqua.res.riseLst1{1, currentEvent}.dlyMap50;
        networkData.propagationMap_event{currentEvent} = propagationMap_event;
    
        % find the cell related to the current event
        cellNumber_event = find(cellfun(@(c) any(c == currentEvent), data_CFU.cfuInfo1(:, 2)));
        networkData.cellNumber_event{currentEvent} = cellNumber_event;

        % find if cell is perivascular or nonPerivascular
        if ismember(cellNumber_event,data_analysis.perivascularCells)
            networkData.cellType_event{currentEvent} = 0;
        else
            networkData.cellType_event{currentEvent} = 2;
        end

        % Get map of the cell of current event
        cellMap_event = data_CFU.cfuInfo1{cellNumber_event, 3};
        networkData.cellMap_event{currentEvent} = cellMap_event;
     
        % Find list of all simultaneous events related to that event and fill up simultaneousMatrix
        simultaneousEvents_current = networkData.simultaneousEvents_all{currentEvent};
      
        for a = 1:length(simultaneousEvents_current)
            simultaneousEvent = simultaneousEvents_current(a);

            % find the propagation matrix map related to the simultaneous event
            propagationMap_all = data_aqua.res.riseLst1{1, simultaneousEvent}.dlyMap50;
            networkData.propagationMap_all{currentEvent}{end+1} = propagationMap_all;
            
            % Find the cell of the specific event
            cellNumber_all = find(cellfun(@(c) any(c == simultaneousEvent), data_CFU.cfuInfo1(:, 2)));

            if ~isempty(cellNumber_all)
                if isempty(networkData.cellNumber_all{currentEvent})
                    networkData.cellNumber_all{currentEvent} = cellNumber_all;
                else
                    networkData.cellNumber_all{currentEvent} = [networkData.cellNumber_all{currentEvent}, cellNumber_all];
                end
            end
            
            % find if cell is perivascular or nonPerivascular
            if ismember(cellNumber_all, data_analysis.perivascularCells)
                cellType_all = 0;
            else
                cellType_all = 2;
            end
            
            % populate networkData table
            if isempty(networkData.cellType_all{currentEvent})
               networkData.cellType_all{currentEvent} = cellType_all;
            else
                networkData.cellType_all{currentEvent} = [networkData.cellType_all{currentEvent}, cellType_all];
            end

            % Find the map of the simultaneousEvent to the cell array
            cellMap_all = data_CFU.cfuInfo1{cellNumber_all, 3};
            networkData.cellMap_all{currentEvent}{end+1} = cellMap_all;  

%             % Plot image of propagation within cell map
%             [networkData, newTempMap] = propagationCellMap(networkData, currentEvent, a);

            % remove 'cols_to_delete' from analysis, it includes eventToDelete and events from multinucleated cells
            if ismember(currentEvent, data_analysis.cols_to_delete) 
                simultaneousMatrixByEvent(currentEvent, :) = 0;
                simultaneousMatrixByEvent(:, currentEvent) = 0;
            elseif ismember(simultaneousEvent, data_analysis.cols_to_delete) 
                simultaneousMatrixByEvent(simultaneousEvent, :) = 0;
                simultaneousMatrixByEvent(:, simultaneousEvent) = 0;
            else
                simultaneousMatrixByEvent(currentEvent, simultaneousEvent) = 1;
                simultaneousMatrixByEvent(currentEvent, currentEvent) = 1; 
            end

            % Calculate delay between current event and 1st simultanous event
            element = simultaneousMatrixByEvent(currentEvent, simultaneousEvent); % Access the current element
            % Perform operations on 'element'
            %fprintf('Element at (%d, %d): %d\n', currentEvent, simultaneousEvent, element);
            startingFrame_currentEvent = startingFrame{currentEvent};
            startingFrame_simultaneousEvent = startingFrame{simultaneousEvent};
            delay = startingFrame_simultaneousEvent-startingFrame_currentEvent;
            if delay < 0
               delay = 0;
            end

            if ismember(currentEvent, data_analysis.cols_to_delete)
                % Nullify rows and columns for currentEvent
                simultaneousMatrixDelayByEvent(currentEvent, :) = 0;
                simultaneousMatrixDelayByEvent(:, currentEvent) = 0;

                for cellToDelete = 1:length(data_CFU.cfuInfo1)
                    if ismember(currentEvent, data_CFU.cfuInfo1{cellToDelete,2}) % Check if targetValue exists in column2{i}
                        cellToDelete_event = data_CFU.cfuInfo1{cellToDelete,1}; % Get the corresponding value from column1
                        break; % Exit loop once found
                    end
                end
            
                % Nullify rows and columns for simultaneousEvent if necessary
            elseif ismember(simultaneousEvent, data_analysis.cols_to_delete)
                    simultaneousMatrixDelayByEvent(simultaneousEvent, :) = 0;
                    simultaneousMatrixDelayByEvent(:, simultaneousEvent) = 0;

                for cellToDelete = 1:length(data_CFU.cfuInfo1)
                    if ismember(simultaneousEvent, data_CFU.cfuInfo1{cellToDelete,2}) % Check if targetValue exists in column2{i}
                        cellToDelete_event = data_CFU.cfuInfo1{cellToDelete,1}; % Get the corresponding value from column1
                        break; % Exit loop once found
                    end
                end
            else
                % Assign delay if neither event is in cols_to_delete
                cellToDelete_event = [];
                simultaneousMatrixDelayByEvent(currentEvent, simultaneousEvent) = delay;
            end
                
                % Check if the cell already contains data
            if isempty(simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all})
                if ~isempty(cellToDelete_event)
                    if cellToDelete_event == cellNumber_event || cellNumber_all
                        simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all} = [];
                    else
                        simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all} = delay;
                    end
                else
                    simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all} = delay;
                end
            else
                if ~isempty(cellToDelete_event)
                    if cellToDelete_event == cellNumber_event || cellNumber_all
                        simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all} = ...
                        [simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all}, 0];
                    else
                        simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all} = ...
                        [simultaneousMatrixDelaybyCell{cellNumber_event, cellNumber_all}, delay];
                    end
                end
            end            
        end

        % Duplicate all events and Remove the current event from all simultaneous events
        simultaneousEvents_network = simultaneousEvents_current; %duplicate 
        simultaneousEvents_network(simultaneousEvents_current == currentEvent) = [];
        % Update the networkData table with the modified list
        networkData.simultaneousEvents_network{currentEvent} = simultaneousEvents_network;

        % Duplicate events
        simultaneousEvents_network_corrected = simultaneousEvents_network;
        networkData.simultaneousEvents_network_corrected{currentEvent} = simultaneousEvents_network_corrected;
          
        % Check if any elements in simultaneousEvents_network_corrected needs to be deleted
        if any(ismember(simultaneousEvents_network_corrected, data_analysis.cols_to_delete))
            % Keep only elements that are not in cols_to_delete
            simultaneousEvents_network_corrected = simultaneousEvents_network_corrected(~ismember(simultaneousEvents_network_corrected, data_analysis.cols_to_delete));
            networkData.simultaneousEvents_network_corrected{currentEvent} = simultaneousEvents_network_corrected;
        end

        % Initialize an array to store indices of simultaneousEvents to remove
        indicesToRemove = [];

        %remove current event from network events
        for b = 1:length(simultaneousEvents_network_corrected)
            simultaneousEvent_network_corrected = simultaneousEvents_network_corrected(b);
    
            % Check if the simultaneous event is in the list of cellEvents
            if ismember(simultaneousEvent_network_corrected, data_CFU.cfuInfo1{cellNumber_event, 2})
                % Add the index to the list of indices to remove
                indicesToRemove = [indicesToRemove, b];
            end

            % find the propagation matrix related to the current event
            propagationMap_network = data_aqua.res.riseLst1{1, simultaneousEvent_network_corrected}.dlyMap50;
            networkData.propagationMap_network{currentEvent}{end+1} = propagationMap_network;

            % Find the cell of the specific event
            cellNumber_network = find(cellfun(@(c) any(c == simultaneousEvent_network_corrected), data_CFU.cfuInfo1(:, 2)));

            if ~isempty(cellNumber_network)
                if isempty(networkData.cellNumber_network{currentEvent})
                    networkData.cellNumber_network{currentEvent} = cellNumber_network;
                else
                    networkData.cellNumber_network{currentEvent} = [networkData.cellNumber_network{currentEvent}, cellNumber_network];
                end
            end

            % find if cell is perivascular or nonPerivascular
            if ismember(cellNumber_network, data_analysis.perivascularCells)
                cellType_network = 0;
            else
                cellType_network = 2;
            end

            % populate networkData table
            if isempty(networkData.cellType_network{currentEvent})
               networkData.cellType_network{currentEvent} = cellType_network;
            else
                networkData.cellType_network{currentEvent} = [networkData.cellType_network{currentEvent}, cellType_network];
            end
            
            % Check if the cellMap entry for simultaneousEvent is empty and initialize if needed
            if isempty(networkData.cellMap_network{currentEvent})
                % Initialize networkData.cellMap{currentEvent} as a cell array of size 1x1
                networkData.cellMap_network{currentEvent} = {};
            end

            % Find the map of the simultaneousEvent to the cell array
            cellMap_network = data_CFU.cfuInfo1{cellNumber_network, 3};
            networkData.cellMap_network{currentEvent}{end+1} = cellMap_network;
        end
        
        % Remove the indices in reverse order to avoid out of range errors
        for i = length(indicesToRemove):-1:1
            networkData.simultaneousEvents_network_corrected{currentEvent}(indicesToRemove(i)) = [];
        end

        % Combine cell maps
        [networkData, combinedMatrix, pairwiseDistances] = combineCellMaps(networkData, currentEvent);
    end

    % Get distribution of distances between pair of cells
    % Extract the relevant column
    shortestDistanceBetweenCenters = networkData.shortestDistanceBetweenCenters;
    
    % Initialize an array to store the averages of the upper triangular values
    upperTriAverages = [];
    upperTriValues = [];

    % Conversion factor from pixels to micrometers (spatial resolution = 0.49 um/pixel = 2.04pixel/um)
    conversionFactor = 0.49;
    
    % Process each matrix
    for i = 1:length(shortestDistanceBetweenCenters)
    currentMatrix = shortestDistanceBetweenCenters{i};

        if ismember(i, data_analysis.cols_to_delete)
            shortestDistanceBetweenCenters{i} = [];
        else 
            if isnumeric(currentMatrix)
                % Convert each element from pixels to micrometers
                convertedMatrix = currentMatrix * conversionFactor;
                
                % Extract the upper triangular part, excluding the diagonal
                upperTriValue = convertedMatrix(triu(true(size(convertedMatrix)), 1));
    
                % Append values with no average
                upperTriValues = [upperTriValues; upperTriValue];
            end
        end
    end   

    % Store the upperTriValues for this experiment
    allUpperTriValues{experiment} = upperTriValues;
        
%     % Plot the distribution of the averages
%     figure;
%     %histogram(upperTriAverages);
%     %histogram(upperTriAvg);
%     histogram(upperTriValues);
%     title('Distance between pair of cells of concurrent events (um)');
%     xlabel('Average Distance (µm)');
%     ylabel('Frequency');
% 
%     % Save figure
%     % Edit the filename to end with '_cellDistances'
%     currentFolder = pwd;
%     subfolderName = 'cellDistances'; % Define the subfolder name
%     subfolderPath = fullfile(currentFolder, subfolderName); % Create the full path for the subfolder
%     fileTemp = extractBefore(AquA_fileName, "_AQuA2"); 
%     distancesFileName = strcat(fileTemp, '_cellDistances', '.fig');
%     save_path = fullfile(subfolderPath, distancesFileName);
%     saveas(gcf, save_path);
% 
%     % Save data
%     networkFilename = strcat(fileTemp, '_network_propagation.mat');
%     save(networkFilename);

    % Get events duration and simultaneous events
    eventDuration = table2cell(data_analysis.resultsRaw("Curve - Duration 10% to 10%","ftsTb"));
    eventDuration = eventDuration{1,1}';
    simultaneousEvents_experiment = networkData.simultaneousEvents_network_corrected(:);
    eventDuration_simultaneousEvents = [eventDuration,simultaneousEvents_experiment]; 
    eventDuration_simultaneousEvents_all{experiment} = eventDuration_simultaneousEvents; % Store values for this experiment
    % remove cols_to_delete from analysis
    if ~isempty(data_analysis.cols_to_delete)
        eventDuration_simultaneousEvents_all{experiment}(data_analysis.cols_to_delete, :) = [];
    end

    % MATRIX for all cells (distance and delay)
    pairwiseDistanceMatrix = nan(numCells, numCells); % Initialize with NaN

    % Calculate distance and average delay between events within cells
    for cell1 = 1:numCells
        for cell2 = 1:numCells

            % Retrieve the cell numbers
            cellNumber_matrix = [cell1,cell2]; 

            % Retrieve the cell maps for the two cells
            cellMap1 = data_CFU.cfuInfo1{cell1, 3};
            cellMap2 = data_CFU.cfuInfo1{cell2, 3};
            cellMap_matrix = {cellMap1,cellMap2};
 
            % Initialize the combined matrix with zeros
            combinedMatrix_allCells = zeros(size(cellMap1));
        
            % Loop over each cellMap and combine non-zero values
            for k = 1:numel(cellMap_matrix)
                currentCellMap = cellMap_matrix{k};
                % Create binary mask for non-zero elements
                currentMask = currentCellMap ~= 0;
                
                % Combine non-zero values into the combinedMatrix
                combinedMatrix_allCells(currentMask) = currentCellMap(currentMask);
            end
%             figure;
%             imshow(combinedMatrix_allCells, []); % Display the combined matrix
%             title('Combined Non-Zero Values from All Cell Maps');

            % Compute pairwise distance 
            distance = computePairwiseCenterDistances({cellMap1, cellMap2});
            % Logical indexing to exclude NaN values
            nonNaNValues = distance(~isnan(distance));
            firstValue = nonNaNValues(1); 
            convertedValue = firstValue * conversionFactor;
            
            % Fill the matrix symmetrically
            pairwiseDistanceMatrix(cell1, cell2) = convertedValue;
            pairwiseDistanceMatrix(cell2, cell1) = convertedValue; % Symmetric
 
            % Delay
            % Get the delay data for the current cell
            data = simultaneousMatrixDelaybyCell{cell1,cell2};
            % Exclude zeros and negative values using logical indexing
            validData = data(data > 0);
            % Calculate the mean of the valid data (non-zero and non-negative elements)
            if ~isempty(validData) % Check if there are valid elements
                meanValue = mean(validData);
            else
                meanValue = NaN; % Assign NaN if no valid data
            end    
            % Store the mean value in the output matrix
            simultaneousMatrixDelaybyCell_average{cell1, cell2} = meanValue;
        end
    end

    % Digraph
    [adjMatrix, centers_allCells, G, rowsWithSingleAppearance, nodeDegree] = plotCellDistanceNetwork(data_CFU, data_analysis, simultaneousMatrixDelaybyCell, simultaneousMatrixDelaybyCell_average, pwd, AquA_fileName);

    % Compute nodeOUTdegree x listOfEvents_perCell
    listOfEvents_perCell = cellfun(@numel, data_CFU.cfuInfo1(:,2));
    listOfEvents_perCell_nodeOUTdegree = [listOfEvents_perCell, nodeDegree];
    listOfEvents_perCell_nodeOUTdegree_all{experiment} = listOfEvents_perCell_nodeOUTdegree;
    
    % Plot the resultant vector
    resultantVector = sumEdgeVectors(G, centers_allCells, rowsWithSingleAppearance, pwd, AquA_fileName);

    % Identify connected networks
    saveConnectedNetworks(G, centers_allCells, pwd, AquA_fileName);

    % Compute cell pair distance x edge/delay
    cellPairs_edges_distanceMicron_multipleAppearance = getEdgeDistances(adjMatrix, centers_allCells, numCells, rowsWithSingleAppearance);
    cellPairs_edges_distanceMicron_multipleAppearance_all = [cellPairs_edges_distanceMicron_multipleAppearance_all; cellPairs_edges_distanceMicron_multipleAppearance];

    % Save network data
    fileTemp = extractBefore(AquA_fileName, "_AQuA2");
    pathTemp = extractBefore(pwd, "3."); 
    subfolderNetworkName = '4._network_propagation.mat'; % Define the subfolder name
    subfolderNetworkPath = fullfile(pathTemp, subfolderNetworkName); % Create the full path for the subfolder
    % Create the full file name with path
    networkFilename = fullfile(subfolderNetworkPath, strcat(fileTemp, '_network_propagation.mat'));
    save(networkFilename, '-v7.3');
end

% Correlation - duration of events x number of simultaneous events
correlationBetweenEventDurationAndNumberSimultaneousEvents = computeEventDurationCorrelation(eventDuration_simultaneousEvents_all);

% Distribution - distances between cells with concurrent activity
cleanDataDistances_um = cellPairsDistanceDistribution(allUpperTriValues);

% Correlation - cell pair distance x edge
correlationBetweenCellPairsAndEdge = correlationBetweenCellPairsAndEdge(cellPairs_edges_distanceMicron_multipleAppearance_all);

% Correlation - list of events x node OUTdegree
correlationBetweenNumberofEventsAndOUTdegreeNode = correlationBetweenNumberofEventsAndOUTdegreeNode(listOfEvents_perCell_nodeOUTdegree_all);

