%% Load .mat files 

clear all;

fullCraniotomyCSDDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick\0.49resolution_correct\3._analysisByEvent.mat';
%fullCraniotomyBIBNDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\BIBN\3._analysisByEvent.mat';

cd (fullCraniotomyCSDDir);

% Get all .mat files in the directory
FilesAll = dir(fullfile(fullCraniotomyCSDDir, '*_analysisByEvent.mat')); 

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

%% Extract and combine resultData from each experiment

combinedTable = [];

% Loop through each file
for i = 1:length(sortedFileNames)
    % Load the .mat file
    data = load(sortedFileNames{i});
    
    % Extract data from the structure (assuming variable names are consistent across files)
    % Here, assuming the variable name is 'result' in each .mat file
    resultData = data.resultsFinal;

    % Create a column with the filename
    fileNameColumn = repelem(sortedFileNames(i), size(resultData, 1), 1); % Repeat filename for each row
  
    % Append resultData to combinedTable
    combinedTable = [combinedTable; resultData, table(fileNameColumn)];
end


%% Extract and combine resultData from each experiment

variablesNames = {'Area', ...
            'Perimeter', ...
            'Circularity', ...
            'Max_dFF', ...
            'Duration50to50', ...
            'Duration10to10', ...
            'RisingDuration10to90', ...
            'DecayingDuration90to10',...
            'datAUC', ...
        	'dFFAUC',...
            'NumberOfEvents'};

paramTables_all_byCell = [];
paramTables_all_sum_byFOV = [];
paramTables_all_mean_byFOV = [];
paramTables_all_median_byFOV = [];
eventsByCell_all = [];
eventsDelays = [];
columnsToMedian_all = [];

% Loop through each file
for file = 1:length(sortedFileNames)

    % Load the .mat file
    data_analysis = load(sortedFileNames{file});

    startingFrame = data_analysis.resultsRaw{3,:};

    for eventFrame = 1:size(startingFrame,2)
        if ismember(eventFrame, data_analysis.cols_to_delete)
            eventDelay = [];
        else
            eventDelay = 1854 - startingFrame{1, eventFrame};
        end
        eventsDelays = [eventsDelays; eventDelay];
    end

    % frames for CSD phases
    baseline_preCSD_frames = [1,927]; %[927,1854] [1,927]
    preCSD_frames = 1854;
    % 1min duringCSD
    duringCSD_frames = [1855, 1917]; %60sec (ca. 62 frames)
    postCSD_frames = [1918, 3772];
    % 2min duringCSD
%     duringCSD_frames = [1855, 1979]; %60sec (ca. 124 frames)
%     postCSD_frames = [1980, 3837];

    % separate events in pre, during and postCSD
    eventsByCell_experiment = eventsCSDphases(data_analysis, preCSD_frames, duringCSD_frames, postCSD_frames, baseline_preCSD_frames);
    eventsByCell_all = [eventsByCell_all; eventsByCell_experiment];

    % get median of parameters by CSD phase
    paramTables = CSDparams(eventsByCell_experiment, data_analysis, variablesNames);
    paramTables_all_byCell = [paramTables_all_byCell; paramTables];

    % sum, mean and median of number of events by FOV
    % sum
    [columnsToSum, eventSums] = calculateEventSum(paramTables);
    paramTables_all_sum_byFOV = [paramTables_all_sum_byFOV; eventSums];
    % mean
    [columnsToMean, eventMeans] = calculateEventMean(paramTables);
    paramTables_all_mean_byFOV = [paramTables_all_mean_byFOV; eventMeans];
    % median
    [columnsToMedian, eventMedians] = calculateEventMedian(paramTables);
    paramTables_all_median_byFOV = [paramTables_all_median_byFOV; eventMedians];

    columnsToMedian_all = [columnsToMedian_all; columnsToMedian];
end

% Calculate event rate byCell and byFOV(median)
[eventHz_byCell, eventHz_byFOV] = calculateEventRate(paramTables_all_byCell, paramTables_all_median_byFOV);

%get name of columns/parameters
parameters = fieldnames(paramTables_all_byCell);

% Initialize the output structures
paramTables_allPhases = struct();
paramTables_allPhases_cleanedData = struct();
paramTables_pre_duringCSD_cleanedData = struct();
paramTables_pre_postCSD_cleanedData = struct();
paramTables_during_postCSD_cleanedData = struct();
classifications = struct();

% Directory to save the CSV files
outputDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick';
%outputDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\BIBN';

for param = 1:length(parameters)
    parameter = parameters{param};
    disp(parameter);
    

    % Initialize an empty table to store concatenated data
    combinedTable = table(); % Initialize an empty table to store concatenated data

    for expt = 1:length(paramTables_all_byCell)
        currentTable = paramTables_all_byCell(expt).(parameter);
        combinedTable = [combinedTable; currentTable];
    end

    for x = size(combinedTable, 1):-1:1  % Looping backward to avoid skipping rows
        if isempty(combinedTable.preCSD{x}) && isempty(combinedTable.duringCSD{x}) && isempty(combinedTable.postCSD{x})
            combinedTable(x, :) = [];  % Remove the entire row
        end
    end

    % Store the concatenated table in the output struct
    paramTables_allPhases.(parameter) = combinedTable;

    % use this combinedTable for prism (preCSD_allFeatures), param-by-param
    combinedTable = table2cell(combinedTable);
    emptyRows = cellfun(@isempty, paramTables_allPhases.(parameter){:,2});
    combinedTable(emptyRows, :) = [];

    % Clean and Classify data into increase, decrease, noChange
    [cleanedData, classifiedData, classifiedTable.(parameter)] = processAndClassify(paramTables_allPhases.(parameter));


    % Store cleaned data
    paramTables_allPhases_cleanedData.(parameter) = cleanedData.allPhases;
    paramTables_pre_duringCSD_cleanedData.(parameter) = cleanedData.pre_duringCSD;
    paramTables_pre_postCSD_cleanedData.(parameter) = cleanedData.pre_postCSD;
    paramTables_during_postCSD_cleanedData.(parameter) = cleanedData.during_postCSD;

    % Store classifications
    classifications.(parameter) = classifiedData;

    %writetable(classifiedTable.(parameter).allPhases, fullfile(outputDir, [parameter '_pre_during_postCSD.csv']));
    %writetable(classifiedTable.(parameter).pre_duringCSD, fullfile(outputDir, [parameter '_pre_duringCSD.csv']));
    %writetable(classifiedTable.(parameter).pre_postCSD, fullfile(outputDir, [parameter '_pre_postCSD.csv']));
    %writetable(classifiedTable.(parameter).during_postCSD, fullfile(outputDir, [parameter '_during_postCSD.csv']));

    %savePath = extractBefore(fullCraniotomyCSDDir, '\3.');
    %plotCSDBoxplots(paramTables_allPhases_cleanedData, parameter, savePath);
end

%% cell count for slide 66

cellType = combinedTable{:, 6}; % Extract column 6 as an array
cellType = cell2mat(cellType);
totalPerivascular = sum(cellType == 0); % Count rows with 0 or 2
totalNonPerivascular = sum(cellType == 2);

cellType_pre_duringCSD = paramTables_pre_duringCSD_cleanedData.Area(:,3); % Extract column 6 as an array
pre_duringCSD_totalPerivascular = sum(cellType_pre_duringCSD == 0); % Count rows with 0 or 2
pre_duringCSD_totalNonPerivascular = sum(cellType_pre_duringCSD == 2);

cellType_pre_postCSD = paramTables_pre_postCSD_cleanedData.Area(:,3); % Extract column 6 as an array
pre_duringCSD_totalPerivascular = sum(cellType_pre_postCSD == 0); % Count rows with 0 or 2
pre_duringCSD_totalNonPerivascular = sum(cellType_pre_postCSD == 2);

combinedTable2 = combinedTable;

for x = size(combinedTable2, 1):-1:1  % Looping backward to avoid skipping rows
    if isempty(combinedTable2.baseline_preCSD{x})
        combinedTable2(x, :) = [];  % Remove the entire row
    end
end

%% to compare with baseline data

baselineCSD_area = paramTables_allPhases.Area.baseline_preCSD(~cellfun(@isempty, paramTables_allPhases.Area.baseline_preCSD));
baselineCSD_perimeter = paramTables_allPhases.Perimeter.baseline_preCSD(~cellfun(@isempty, paramTables_allPhases.Perimeter.baseline_preCSD));
baselineCSD_circularity = paramTables_allPhases.Circularity.baseline_preCSD(~cellfun(@isempty, paramTables_allPhases.Circularity.baseline_preCSD));
baselineCSD_maxDFF = paramTables_allPhases.Max_dFF.baseline_preCSD(~cellfun(@isempty, paramTables_allPhases.Max_dFF.baseline_preCSD));
baselineCSD_duration10to10 = paramTables_allPhases.Duration10to10.baseline_preCSD(~cellfun(@isempty, paramTables_allPhases.Duration10to10.baseline_preCSD));
baselineCSD_AUCdFF = paramTables_allPhases.dFFAUC.baseline_preCSD(~cellfun(@isempty, paramTables_allPhases.dFFAUC.baseline_preCSD));

% baseline_preCSD (15min)
baseline_preCSD_eventList = eventsByCell_all(:,5);
baseline_preCSD_eventList_cleaned = baseline_preCSD_eventList(~cellfun(@isempty, baseline_preCSD_eventList));
baseline_preCSD_numberOfEvents = cellfun(@(x) numel(x), baseline_preCSD_eventList_cleaned);
baseline_preCSD_numberOfEvents_Hz = baseline_preCSD_numberOfEvents/900;
baseline_preCSD_numberOfEvents = [baseline_preCSD_numberOfEvents, baseline_preCSD_numberOfEvents_Hz];

% baseline overall (15min)
baseline_numberOfEvents = [numberOfEvents_perivascular_NM; numberOfEvents_nonPerivascular_NM];
baseline_numberOfEvents = table2cell(baseline_numberOfEvents);
baseline_numberOfEvents = cell2mat(baseline_numberOfEvents);
baseline_numberOfEvents_Hz = baseline_numberOfEvents/900;

%% number of events by phase

% preCSD
preCSD_eventList = eventsByCell_all(:,2);
preCSD_eventList_cleaned = preCSD_eventList(~cellfun(@isempty, preCSD_eventList));
preCSD_numberOfEvents = cellfun(@(x) numel(x), preCSD_eventList_cleaned);
preCSD_numberOfEvents_Hz = preCSD_numberOfEvents/1800;
preCSD_numberOfEvents = [preCSD_numberOfEvents, preCSD_numberOfEvents_Hz];

% duringCSD
duringCSD_eventList = eventsByCell_all(:,3);
duringCSD_eventList_cleaned = duringCSD_eventList(~cellfun(@isempty, duringCSD_eventList));
duringCSD_numberOfEvents = cellfun(@(x) numel(x), duringCSD_eventList_cleaned);
%duringCSD_numberOfEvents_Hz = duringCSD_numberOfEvents/1800;
%duringCSD_numberOfEvents = [duringCSD_numberOfEvents, duringCSD_numberOfEvents_Hz];

% postCSD
postCSD_eventList = eventsByCell_all(:,4);
postCSD_eventList_cleaned = postCSD_eventList(~cellfun(@isempty, postCSD_eventList));
postCSD_numberOfEvents = cellfun(@(x) numel(x), postCSD_eventList_cleaned);
postCSD_numberOfEvents_Hz = postCSD_numberOfEvents/1800;
postCSD_numberOfEvents = [postCSD_numberOfEvents, postCSD_numberOfEvents_Hz];

%% stackedBar for number of events 

% Data for the categories (based on manual count of responses on each FOV: eventHz_byCell)
stackedBarCSD = [
    19, 19, 16, 8, 15, 3, 4, 16, 17, 18, 11, 7, 11;  % Decrease
    1, 0, 14, 4, 4, 2, 3, 2, 1, 1, 5, 7, 5;  % Increase
    1, 2, 4, 5, 3, 1, 6, 9, 0, 4, 0, 0, 1  % No Change
];

stackedBarBIBN = [
    14, 3, 14;  % Decrease
    3, 4, 3;  % Increase
    1, 0, 0  % No Change
];

% Categories
categories = {'Decrease', 'Increase', 'No Change'};

% Create stacked bar chart
figure;
bar(stackedBarBIBN', 'stacked');

% Add labels and title
xlabel('FOV');
ylabel('# cells');
title('Cellular Event Rate (Hz) Response: DuringCSD vs. PreCSD');
legend(categories, 'Location', 'best');

% Set maximum value for y-axis
ylim([0 20]);  % Adjust the maximum value (e.g., 20) as needed

grid off;

%% pie charts

% Data for the categories (rows are categories, columns are FOVs)
stackedBarCSD = [
    19, 19, 16, 8, 15, 3, 4, 16, 17, 18, 11, 7, 11;  % Decrease
    1, 0, 14, 4, 4, 2, 3, 2, 1, 1, 5, 7, 5;  % Increase
    1, 2, 4, 5, 3, 1, 6, 9, 0, 4, 0, 0, 1  % No Change
];

stackedBarBIBN = [
    14, 3, 14;  % Decrease
    3, 4, 3;  % Increase
    1, 0, 0  % No Change
];

% Normalize each category by its total to calculate percentages
percentIncrease = stackedBarBIBN(1, :) / sum(stackedBarBIBN(1, :)) * 100;
percentDecrease = stackedBarBIBN(2, :) / sum(stackedBarBIBN(2, :)) * 100;
percentNoChange = stackedBarBIBN(3, :) / sum(stackedBarBIBN(3, :)) * 100;

% Sort each category in descending order
[percentIncreaseSorted, idxIncrease] = sort(percentIncrease, 'descend');
[percentDecreaseSorted, idxDecrease] = sort(percentDecrease, 'descend');
[percentNoChangeSorted, idxNoChange] = sort(percentNoChange, 'descend');

% Generate corresponding labels
labelsIncrease = arrayfun(@(x) ['FOV ' num2str(x)], idxIncrease, 'UniformOutput', false);
labelsDecrease = arrayfun(@(x) ['FOV ' num2str(x)], idxDecrease, 'UniformOutput', false);
labelsNoChange = arrayfun(@(x) ['FOV ' num2str(x)], idxNoChange, 'UniformOutput', false);

% Create pie charts
figure;

% Pie chart for "Increase"
subplot(1, 3, 1);
pie(percentIncreaseSorted, labelsIncrease);
title('Increase Event Rate (Hz)', 'FontSize',15);
%legend(arrayfun(@(x) ['FOV ' num2str(x)], 1:size(stackedBar, 2), 'UniformOutput', false), 'Location', 'bestoutside');

% Pie chart for "Decrease"
subplot(1, 3, 2);
pie(percentDecreaseSorted, labelsDecrease);
title('Decrease Event Rate (Hz)', 'FontSize',15);
%legend(arrayfun(@(x) ['FOV ' num2str(x)], 1:size(stackedBar, 2), 'UniformOutput', false), 'Location', 'bestoutside');

% Pie chart for "No Change"
subplot(1, 3, 3);
pie(percentNoChangeSorted, labelsNoChange);
title('No Change Event Rate (Hz)', 'FontSize',15);
%legend(arrayfun(@(x) ['FOV ' num2str(x)], 1:size(stackedBar, 2), 'UniformOutput', false), 'Location', 'bestoutside');