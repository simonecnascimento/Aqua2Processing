%% AQuA2_CFU.m read metadata file from AQuA2 and CFU GUI. 
% For all events detected, it gives 15 values/features, including: Index number (=event number), in sequential order,
% area, circularity, max Dff, duration, etc, for all events. The new AQuA2 also outputs information about detected cells ("CFU") which 
% is provided as cfuInfo1/2, depending on the channel chosen. It gives the list of cells detected and their corresponding list of events.
% This detection has to be somewhat manual, in order to adjust overlap threshold (important!!!) and identify areas where detection is not clear

clear all;

spreadsheetPath = 'R:\Levy Lab\2photon\ImagingDatasets_Simone_231218.xlsx';
sheet = 'Macrophage';
dataTable = readcell(spreadsheetPath, 'sheet',sheet);  
colNames = dataTable(1,:); %remove all info and leave just column names
dataTable(1,:) = []; %remove column names

dataCol = struct( ...
    'mouse',find(contains(colNames, 'Mouse')), ...
    'date',find(contains(colNames, 'Date')), ...
    'FOV',find(contains(colNames, 'FOV')), ...
    'run',find(contains(colNames, 'Runs')), ...
    'Zproj',find(contains(colNames, 'Zproj')));

% Convert the numeric values in the date column to strings, regardless of their original data type
dataTable(:,dataCol.date) = cellfun(@num2str, dataTable(:,dataCol.date), 'UniformOutput',false); 

% --- Specify the row number of interest
rowOfInterest = 24; 

% Extract information from the specified row
mouse = dataTable{rowOfInterest, dataCol.mouse};
date = dataTable{rowOfInterest, dataCol.date};
fov = dataTable{rowOfInterest, dataCol.FOV};
run = dataTable{rowOfInterest, dataCol.run};
zProj = dataTable{rowOfInterest, dataCol.Zproj};

%% Load metadata files generated by AQuA

folderPath = fullfile('D:\2photon\Simone\Simone_Macrophages\', mouse, '\', [date '_FOV' num2str(fov)], '\AQuA2\1Hz');
%substackOptions = {'Substack(1-900)', 'Substack(901-1800)'};  % Add more options if needed

% List the contents of the folder
contents = dir(folderPath);

% Exclude '.', '..', and 'Thumbs.db' entries, so 'contents' contains only the actual files and folders in the directory
excludeList = {'.', '..', 'Thumbs.db'};
contents = contents(~ismember({contents.name}, excludeList));

% Initialize subfolder name
subfolderName = '';

% Loop through the contents and find the subfolder with 'Substack' in its name
for s = 1:numel(contents)
    if contents(s).isdir && contains(contents(s).name, 'Substack') 
        subfolderName = contents(s).name;
        disp(subfolderName);

        % Construct the file path for CFU based on the information from the spreadsheet
        AQuA2_CFUinfo_filePath = fullfile( ...
        'D:\2photon\Simone\Simone_Macrophages\', mouse, '\', [date '_FOV' num2str(fov)], '\AQuA2\1Hz', [subfolderName, '_AQuA_res_cfu']);
        AQuA2_CFUinfo = load(AQuA2_CFUinfo_filePath);

        % Construct the file path for AQuA results based on the information from the spreadsheet
        AQuA2_Feature_filePath = fullfile( ...
        'D:\2photon\Simone\Simone_Macrophages\', mouse, '\', [date '_FOV' num2str(fov)], '\AQuA2\1Hz', ...
        [subfolderName, '\'], [subfolderName, '_AQuA2']);
        AQuA2_Feature = load(AQuA2_Feature_filePath);

        % Create a table (finalResults) with final results of all cells and median of features
        %Features to remove 
        featuresIndicesToDelete = [1,3,7,15:24];
        
        %Get metadata results of ALL events
        resultsRaw = AQuA2_Feature.res.featureTable1;
        resultsExpanded = cell2table(resultsRaw.ftsTb); %expanded events without features column
        resultsExpanded(featuresIndicesToDelete,:) = [];
        
        %Create a resultsFinal table of features
        resultsFinal = resultsRaw; %copy of metadata results table
        resultsFinal.ftsTb = []; %assign empty values
        resultsFinal(featuresIndicesToDelete,:) = []; %remove features

        %Add new features to resultsFinal table    
        existingRowNames = resultsFinal.Properties.RowNames; % Getting the existing row names
        newRowNames = {'Number of Events'; 'Cell location'; 'Red label'; 'dFF'};
        updatedRowNames = [existingRowNames; newRowNames];
        newRowsData = zeros(4, width(resultsFinal));   % Adding new rows to the table (e.g., with zeros)
        resultsFinal = [resultsFinal; array2table(newRowsData, 'VariableNames', resultsFinal.Properties.VariableNames)];
        resultsFinal.Properties.RowNames = updatedRowNames; % Setting the updated row names to the table
        
        % Results of CFU
        % Initialize variables and create a table
        cellList = AQuA2_CFUinfo.cfuInfo1(:,1); 
        eventList = AQuA2_CFUinfo.cfuInfo1(:,2);
        cellList_eventList = table(cellList, eventList);
        
        % Transform table to struct array
        cellList_eventList = table2struct(cellList_eventList);
        
        % Ask the user to enter multiple sets of cells to merge
        cellsToMerge = input('Enter ID of cells to merge (e.g., [1,2; 3,4]): ');
        
        % Ask the user to enter ID of perivascular cells
        perivascularCells = input('Enter ID of perivascular cells (e.g., [1,2,3,4]): ');

        % Ask the user to enter ID of cells adjacent to the vasculature
        adjacentCells = input('Enter ID of cells adjacent to the vasculature (e.g., [1,2,3,4]): ');

        % Ask the user to enter ID of cells labelled with red
        redLabel = input('Enter ID of cells labelled with red (e.g., [1,2,3,4]): ');

        % Check if userInput is valid (you may want to add more validation)
        if isnumeric(cellsToMerge) && ismatrix(cellsToMerge)                          
        % Iterate through cellsToMerge
            for m = 1:size(cellsToMerge, 1)
                % Extract cell indices to merge
                cellIndices = cellsToMerge(m, :);
        
                % Combine events of CFUs
                combinedCellEvents = cat(1, cellList_eventList(cellIndices).eventList);
            
                % Update the first CFU with the combined events
                cellList_eventList(cellIndices(1)).eventList = combinedCellEvents;
            end

            % Remove the extra CFUs
            cellList_eventList(cellsToMerge(:,2:end)) = [];            
                   
            % End of the block for processing user input for each set of cells
        else
            disp('Invalid input. Please enter valid cell indices.');
            cellsToMerge = input('Enter pairs of cells to merge (e.g., [1,2; 3,4]): ');
        end             
               
        % Ask the user to enter the events to remove from the analysis
        eventsToDelete = input('Enter events to remove (e.g., [1,2,3,4]):');

        % Iterate through CFU cell array
        for c = 1:size(cellList_eventList,1)
            
            % From one cell, extract its ID and its list of events - do for all cells
            cellID = cellList_eventList(c).cellList;
            cellEvents = cellList_eventList(c).eventList;
            numberOfEvents = size(cellEvents,1);
           
            % Initialize an array to store median values for a group of events          
            EventFeatures = table2array(resultsExpanded(:,cellEvents));        
           
            if ~isempty(eventsToDelete)
                for j = 1:numel(eventsToDelete)
                    eventToDelete = eventsToDelete(j);                              
    
                    % Check if the event exists in cellEvents{1, 1}
                    if ismember(eventToDelete, cellEvents)
                        % Find and remove columns in EventFeatures corresponding to events to delete
                        columnsToDelete = ismember(cellEvents, eventToDelete);
                        EventFeatures(:, columnsToDelete) = [];
                    else
                        disp(['Event ', num2str(eventToDelete), ' not found for Cell ', num2str(c)]);
                    end
                end
            end
            
            %Calculate median of the remaining events
            EventFeatures_median = array2table(median(EventFeatures,2),'VariableNames',{['cell ' num2str(cellID) '_median']});

            % Pad EventFeatures_median with NaNs to make it have the same number of rows as resultsFinal
            missingRows = size(resultsFinal, 1) - size(EventFeatures_median, 1);
            emptyRows = array2table(NaN(missingRows, width(EventFeatures_median)), 'VariableNames', EventFeatures_median.Properties.VariableNames);
            EventFeatures_median = [EventFeatures_median; emptyRows];
            
            % Populate resultFinal table with median of rows 1-11
            resultsFinal = [resultsFinal, EventFeatures_median]; %median of each cell, all features, iterating cell in cellList
            
            % Populate resultFinal table with new rows 

            %Row 12 - number of events per cell
            resultsFinal{12,c} = numberOfEvents;

            %Row 13 - cell location
            if c == perivascularCells
                resultsFinal{13,c} = 1;
            elseif c == adjacentCells
                resultsFinal{13,c} = 2;
            else
                resultsFinal{13,c} = 0;
            end 

            %Row 14 - labelled with red dextran
            if c == redLabel
                resultsFinal{14,c} = 1;
            else
                resultsFinal{14,c} = 0;
            end 

            %Row 15 - dFF curve
            dFFcurve = 
                resultsFinal{14,c} = 1;
            else
                resultsFinal{14,c} = 0;
            end 


%           %Average curve of CFU ---- remove events before  ----- ISSUE         
%           cfuCurve = AQuA2_CFUinfo.cfuInfo1{c,5};
%           resultsFinal(13,c) = cfuCurve; 
% 
%           %Plot the average curve of all CFUs
%           CFU_curve = figure;
%           plot(cfuInfo1{3,5})

            %Combine curves
            %Save figure
        end

        %Adjustments to Results table before saving
        resultsFinal = resultsFinal(2:end,:); %remove Index row (Event number)               
        resultsFinal = resultsFinal(:, any(~isnan(table2array(resultsFinal)), 1)); %remove column in which cell has 0 events (e.g: event was remove but cell was not merged. ie. cell was an artifact)
                     
        % Save results file for single experiment
        % Set the output directory to save files
        outputDirFOV = fileparts(AQuA2_CFUinfo_filePath);
        cd(outputDirFOV); %change directory of current folder                
        outputDirAnimal = fullfile('D:\2photon\Simone\Simone_Macrophages\', mouse, '\', 'AQuA2_Results\'); %set directory of animal folder for AQuA2 results

        [filepath,name,ext] = fileparts(AQuA2_CFUinfo_filePath);
        fileTemp = extractBefore(name, "_AQuA_res_cfu");
        save(fullfile(outputDirFOV, strcat(fileTemp,'_analysis'))); % save metadata inside FOV folder       
        writetable(resultsFinal,fullfile(outputDirFOV,strcat(fileTemp,"_resultsFinal")),"WriteRowNames",true,"WriteVariableNames", true, "FileType","spreadsheet"); %save resultsFinal inside FOV folder
        writetable(resultsFinal,fullfile(outputDirAnimal,strcat(fileTemp,"_resultsFinal")),"WriteRowNames",true,"WriteVariableNames", true, "FileType","spreadsheet"); %save resultsFinal to animal folder
    end
end

%% CFU map
%Show all cells in a grid figure - Assuming your cell array is named 'cfuInfo1'

% Determine the number of cells in the array
numCells = size(cfuInfo1, 1);

% Determine the number of rows and columns for the subplot arrangement
numRows = ceil(sqrt(numCells));
numColumns = ceil(numCells / numRows);

% Create a new figure for all subplots
figure;

for cellIndex = 1:numCells
    % Access the cell map based on the current iteration
    cellMap = cfuInfo1{cellIndex, 3};

    % Create a logical mask for non-zero values
    nonZeroMask = cellMap > 0;

    % Create subplots in a grid
    subplot(numRows, numColumns, cellIndex);

    % Display the cell map with only non-zero values
    imshow(nonZeroMask);

    % Add a title to each subplot
    title(['Cell ', num2str(cellIndex)]);
end

% Save the figure
saveas(gcf, 'combined_cells_figure.png');
    