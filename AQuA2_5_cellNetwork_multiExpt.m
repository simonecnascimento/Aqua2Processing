%% Load .mat files  (FULL CRANIOTOMY BASELINE)

clear all;

% Change you Current Folder - fullCraniotomy
cd D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\4._analysisByEvent.mat;

% Initialize an empty table
table = table();

% fullCraniotomy
FilesAll = {
'Pf4Ai162-2_221130_FOV6_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-10_230628_FOV1_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-10_230628_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-10_230628_FOV3_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',... 
'Pf4Ai162-10_230628_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-10_230628_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',... 
'Pf4Ai162-10_230628_FOV7_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-10_230628_FOV8_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV1_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV3_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV6_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV8_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-11_230630_FOV9_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV1_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV3_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV7_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV8_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV9_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-12_230717_FOV10_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-13_230719_FOV1_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-13_230719_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-13_230719_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-13_230719_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-13_230719_FOV6_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-13_230719_FOV7_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-13_230719_FOV8_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-18_240221_FOV1_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-20_240229_FOV1_reg_green_Substack(1-927)_analysisByEvent.mat'};

%% Load .mat files  (THIN BONE BASELINE)

% Change you Current Folder - thinBone
cd D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\thinBone\_analysisByEvent.mat

% Initialize an empty table
table = table();

% thinBone
FilesAll = {
'Pf4Ai162-5_230502_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-6_230502_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-6_230502_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV1_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV3_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV7_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-9_230614_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-9_230614_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat'};

%% Load .mat files  (FULL CRANIOTOMY CSD)

% Change you Current Folder - thinBone
cd D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick\1._analysis.mat;

% Initialize an empty table
table = table();

% thinBone
FilesAll = {
'Pf4Ai162-5_230502_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-6_230502_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-6_230502_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV1_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV2_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV3_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-8_230614_FOV7_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-9_230614_FOV4_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat',...
'Pf4Ai162-9_230614_FOV5_run1_reg_Z01_green_Substack(1-927)_analysisByEvent.mat'};


%% extract network spatial density

% Initialize a cell array to store the upperTriAverages for each experiment
allUpperTriValues = cell(length(FilesAll), 1);

for experiment = 1:length(FilesAll)
    
    % Load analysis .mat file
    data_analysis = load(FilesAll{experiment});
    features = data_analysis.resultsRaw.Row;
    results_complete = data_analysis.resultsRaw.Variables;

    % AQuA2 directory 
    fileTemp_parts = strsplit(data_analysis.filename, '_');
    aqua_directory = fullfile('D:\2photon\Simone\Simone_Macrophages\', ...
        fileTemp_parts{1,1}, '\', ...
        [fileTemp_parts{1,2} '_' fileTemp_parts{1,3}], '\AQuA2\', ...
        [fileTemp_parts{1,1} '_' fileTemp_parts{1,2} '_' fileTemp_parts{1,3} '_run1_reg_Z01_green_Substack(1-927)']);
    AquA_fileName = [data_analysis.filename '_AQuA2.mat'];

%     fileTemp_parts = strsplit(data_analysis.filename, '_');
%     aqua_directory = fullfile('D:\2photon\Simone\Simone_Macrophages\', ...
%         fileTemp_parts{1,1}, '\', ...
%         [fileTemp_parts{1,2} '_' fileTemp_parts{1,3}], '\AQuA2\', ...
%         [fileTemp_parts{1,1} '_' fileTemp_parts{1,2} '_' fileTemp_parts{1,3} '_reg_green_Substack(1-927)']);
%     AquA_fileName = [data_analysis.filename '_AQuA2.mat'];

    % Load AQuA2.mat file 
    fullFilePath = fullfile(aqua_directory, AquA_fileName);
    data_aqua = load(fullFilePath);

    % FILL TABLE WITH CELL NUMBERS, MAPS, AND EVENTS PROPAGATIONS

    % get network data from aqua
    nSimultanousEvents = num2cell(data_aqua.res.fts1.networkAll.nOccurSameTime);
    simultaneousEvents_all = data_aqua.res.fts1.networkAll.occurSameTimeList;

    % Define the number of events
    numEvents = height(simultaneousEvents_all);

    % for current event
    propagationMap_event = cell(numEvents,1);
    % for cell corresponding the current event
    cellNumber = cell(numEvents,1);
    cellMap = cell(numEvents,1);
    
    % for all simultaneous events
    propagationMap_all = cell(numEvents,1);
    % for cells corresponding the all simultaneous events
    cellNumber_all = cell(numEvents,1);
    cellMap_all = cell(numEvents,1);
   
    % for events removing current Event
    simultaneousEvents_network = simultaneousEvents_all; %duplicate 
    propagationMap_network = cell(numEvents,1);
    % for cell corresponding to events in the network
    cellNumber_network = cell(numEvents,1);
    cellMap_network = cell(numEvents,1);

    % combine propagation map of event to map of cell
    eventNetwork_cellMap_all = cell(numEvents,1);
    % calculate shortest distance between centers
    shortestDistanceBetweenCenters = cell(numEvents,1);

    % create a table with network results
    networkData = [propagationMap_event, cellNumber, cellMap, nSimultanousEvents, simultaneousEvents_all, propagationMap_all, cellNumber_all, cellMap_all, eventNetwork_cellMap_all, shortestDistanceBetweenCenters, simultaneousEvents_network, propagationMap_network, cellNumber_network, cellMap_network];
    networkData = cell2table(networkData);
    networkData.Properties.VariableNames = {'propagationMap_event','cellNumber', 'cellMap', 'nSimultanousEvents', 'simultaneousEvents_all', 'propagationMap_all', 'cellNumber_all', 'cellMap_all', 'eventNetwork_cellMap_all', 'shortestDistanceBetweenCenters', 'simultaneousEvents_network', 'propagationMap_network', 'cellNumber_network', 'cellMap_network'};
      
%     % Duplicate table to remove repetitive cells in the future
%     networkData_cleaned = networkData;
% %     cellNumbers = num2cell(networkData_cleaned{currentEvent, 'cellNumber_all'}{1}); % Assuming cellNumber_all is a cell array
% %     cellMaps = networkData_cleaned{currentEvent, 'cellMap_all'}{1}; % cellMap_all might already be a cell array

    % CFU directory
    CFU_directory = fullfile('D:\2photon\Simone\Simone_Macrophages\', ...
        fileTemp_parts{1,1}, '\', ...
        [fileTemp_parts{1,2} '_' fileTemp_parts{1,3}], '\AQuA2\');
    CFU_fileName = [data_analysis.filename '_AQuA_res_cfu.mat'];
    
    % Load CFU.mat file 
    CFU_FilePath = fullfile(CFU_directory, CFU_fileName);
    data_CFU = load(CFU_FilePath);

    for currentEvent = 1:size(networkData,1)

        % find the propagation matrix related to the current event
        propagationMap = data_aqua.res.riseLst1{1, currentEvent}.dlyMap50;
        networkData.propagationMap_event{currentEvent} = propagationMap;
    
        % find the cell related to the current event
        cellIndex = find(cellfun(@(c) any(c == currentEvent), data_CFU.cfuInfo1(:, 2)));
        networkData.cellNumber{currentEvent} = cellIndex;
 
        % Get map of the cell of current event
        cellMap = data_CFU.cfuInfo1{cellIndex, 3};
        networkData.cellMap{currentEvent} = cellMap;

        % Find list of all simultaneous events related to that event
        simultaneousEvents_all = networkData.simultaneousEvents_all{currentEvent};
        
        for a = 1:length(simultaneousEvents_all)
            simultaneousEvent_all = simultaneousEvents_all(a);

            % find the propagation matrix related to the simultaneous event
            propagationMap_all = data_aqua.res.riseLst1{1, simultaneousEvent_all}.dlyMap50;
            networkData.propagationMap_all{currentEvent}{end+1} = propagationMap_all;
            
            % Find the cell of the specific event
            cellIndex_all = find(cellfun(@(c) any(c == simultaneousEvent_all), data_CFU.cfuInfo1(:, 2)));

            if ~isempty(cellIndex_all)
                if isempty(networkData.cellNumber_all{currentEvent})
                    networkData.cellNumber_all{currentEvent} = cellIndex_all;
                else
                    networkData.cellNumber_all{currentEvent} = [networkData.cellNumber_all{currentEvent}, cellIndex_all];
                end
            end

            % Find the map of the simultaneousEvent to the cell array
            cellMap_all = data_CFU.cfuInfo1{cellIndex_all, 3};
            networkData.cellMap_all{currentEvent}{end+1} = cellMap_all;  

            % eventNetwork_cellMap_all
   
%             % Cell Map
%             tempMap = networkData.cellMap_all{currentEvent}{simultaneousEvent_all};
%             nonZeroMask = tempMap > 0;
%             countTempMap = nnz((tempMap));
%             imshow(nonZeroMask);
%             
%             % Event Propagation
%             tempPropagation = networkData.propagationMap_all{currentEvent}{simultaneousEvent_all};  
%             countTempPropagation = nnz(~isnan(tempPropagation));
%             
%             % Get the linear indices of non-zero elements in temp
%             nonZeroIndicesTempMap = find(tempMap ~= 0);
%             
%             % Get the linear indices of non-NaN elements in matrix
%             nonNaNIndicesTempPropagation = find(~isnan(tempPropagation));
%             
%             % Ensure both have the same number of elements
%             if length(nonZeroIndicesTempMap) ~= length(nonNaNIndicesTempPropagation)
%                 error('The number of non-zero elements in temp must equal the number of non-NaN elements in matrix.');
%             end
%             
%             % Create a new variable to hold the result
%             newTempMap = tempMap;
%             
%             % Insert the values from matrix into the newTemp at the corresponding positions
%             newTempMap(nonZeroIndicesTempMap) = tempPropagation(nonNaNIndicesTempPropagation);
%             
%             
%             % Replace all zeros with NaN in temp
%             newTempMap(newTempMap == 0) = NaN;
% 
%             % Substitute networkData table
%             networkData.eventNetwork_cellMap_all{currentEvent}{end+1} = newTempMap; 
            
%             % Plot image of eventNetwork_cellMap_all
%             imagesc(newTempMap)
%             matrixTemp = newTempMap;
%             
%             % Find the maximum value
%             maxValue = max(matrixTemp(:));
%             % Find all positions of the maximum value
%             [maxRows, maxCols] = find(matrixTemp == maxValue);
%             
%             % Find the minimum value
%             minValue = min(matrixTemp(:));
%             % Find all positions of the minimum value
%             [minRows, minCols] = find(matrixTemp == minValue);
%             
%             % Display the results for maximum value
%             disp(['Maximum value: ', num2str(maxValue)]);
%             disp('Positions of max value:');
%             for j = 1:length(maxRows)
%                 disp(['(', num2str(maxRows(j)), ', ', num2str(maxCols(j)), ')']);
%             end
%             
%             % Display the results for minimum value
%             disp(['Minimum value: ', num2str(minValue)]);
%             disp('Positions of min value:');
%             for j = 1:length(minRows)
%                 disp(['(', num2str(minRows(j)), ', ', num2str(minCols(j)), ')']);
%             end
        end

        % Duplicate all events and Remove the current event from all simultaneous events
        simultaneousEvents_network = simultaneousEvents_all; %duplicate 
        simultaneousEvents_network(simultaneousEvents_all == currentEvent) = [];
        % Update the networkData table with the modified list
        networkData.simultaneousEvents_network{currentEvent} = simultaneousEvents_network;
          
        % Initialize an array to store indices of simultaneousEvents to remove
        indicesToRemove = [];

        %remove current event from network events
        for b = 1:length(simultaneousEvents_network)
            simultaneousEvent_network = simultaneousEvents_network(b);
    
            % Check if the simultaneous event is in the list of cellEvents
            if ismember(simultaneousEvent_network, data_CFU.cfuInfo1{cellIndex, 2})
                % Add the index to the list of indices to remove
                indicesToRemove = [indicesToRemove, b];
            end

            % find the propagation matrix related to the current event
            propagationMap_network = data_aqua.res.riseLst1{1, simultaneousEvent_network}.dlyMap50;
            networkData.propagationMap_network{currentEvent}{end+1} = propagationMap_network;

            % Find the cell of the specific event
            cellIndex_network = find(cellfun(@(c) any(c == simultaneousEvent_network), data_CFU.cfuInfo1(:, 2)));

            if ~isempty(cellIndex_network)
                if isempty(networkData.cellNumber_network{currentEvent})
                    networkData.cellNumber_network{currentEvent} = cellIndex_network;
                else
                    networkData.cellNumber_network{currentEvent} = [networkData.cellNumber_network{currentEvent}, cellIndex_network];
                end
            end

            % Check if the cellMap entry for simultaneousEvent is empty and initialize if needed
            if isempty(networkData.cellMap_network{currentEvent})
                % Initialize networkData.cellMap{currentEvent} as a cell array of size 1x1
                networkData.cellMap_network{currentEvent} = {};
            end

            % Find the map of the simultaneousEvent to the cell array
            cellMap_network = data_CFU.cfuInfo1{cellIndex_network, 3};
            networkData.cellMap_network{currentEvent}{end+1} = cellMap_network;
        end
        
        % Remove the indices in reverse order to avoid out of range errors
        for i = length(indicesToRemove):-1:1
            networkData.simultaneousEvents_network{currentEvent}(indicesToRemove(i)) = [];
        end

        % COMBINE CELL MAPS
    
        cellNumbers = num2cell(networkData{currentEvent, 'cellNumber_all'}{1}); % Assuming cellNumber_all is a cell array
        cellMaps = networkData{currentEvent, 'cellMap_all'}{1}; % cellMap_all might already be a cell array
    
        if iscell(cellNumbers)
            % Flatten the cell array and find unique numbers
            cellNumbers = [cellNumbers{:}];
        end
        
        % Remove duplicates in cellNumbers
        [uniqueCellNumbers, ia] = unique(cellNumbers, 'stable');
        
        % Extract corresponding cellMaps for unique cellNumbers
        cleanedCellMaps = cellMaps(ia);
    
        % Update the cleaned table
        networkData{currentEvent, 'cellNumber_all'} = {uniqueCellNumbers};
        networkData{currentEvent, 'cellMap_all'}{:} = cleanedCellMaps;
    
        % Initialize the combined matrix with zeros
        combinedMatrix = zeros(size(cleanedCellMaps{1}));
    
        % Loop over each cellMap and combine non-zero values
        for k = 1:numel(cleanedCellMaps)
            currentCellMap = cleanedCellMaps{k};
            % Create binary mask for non-zero elements
            currentMask = currentCellMap ~= 0;
            
            % Combine non-zero values into the combinedMatrix
            combinedMatrix(currentMask) = currentCellMap(currentMask);
        end
%         figure;
%         imshow(combinedMatrix, []); % Display the combined matrix
%         title('Combined Non-Zero Values from All Cell Maps');
    
        networkData.eventNetwork_cellMap_all{currentEvent} = combinedMatrix;
    
        % Compute the shortest distance
        pairwiseDistances = computePairwiseCenterDistances(cleanedCellMaps);
        networkData.shortestDistanceBetweenCenters{currentEvent} = pairwiseDistances;
        % Create heatmap of the distance matrix
        %createHeatmap(pairwiseDistances);
   
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
    
        if isnumeric(currentMatrix)
            % Convert each element from pixels to micrometers
            convertedMatrix = currentMatrix * conversionFactor;
            
            % Extract the upper triangular part, excluding the diagonal
            upperTriValue = convertedMatrix(triu(true(size(convertedMatrix)), 1));
            
%             % Compute the mean of the upper triangular values, ignoring NaNs
%             upperTriAvg = mean(upperTriValue, 'omitnan');
%             % Append the computed average to the array
%             upperTriAverages = [upperTriAverages; upperTriAvg];

            % Append values with no average
            upperTriValues = [upperTriValues; upperTriValue];
        end
    end   
        
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
% 
%     % Store the upperTriAverages for this experiment
%     allUpperTriAverages{experiment} = upperTriAverages;

    % Store the upperTriValues for this experiment
    allUpperTriValues{experiment} = upperTriValues;


end

%% Get all averages of distances and get distribution




% Initialize an empty array to store all values
allValuesDistances_um = [];

% Loop through each array
for k = 1:length(allUpperTriValues) %allUpperTriAverages
    currentVector = allUpperTriValues{k}; %allUpperTriAverages
    
    % Append the values to the allValues array
    allValuesDistances_um = [allValuesDistances_um; currentVector];
   
end

% Remove NaN values
cleanDataDistances_um = allValuesDistances_um(~isnan(allValuesDistances_um));

% Plot the histogram of the cleaned data
figure;
histogram(cleanDataDistances_um);
title('Distance Metrics for Cells with Concurrent Activity');
xlabel('Distance (um)');
ylabel('Number of cell pairs');

