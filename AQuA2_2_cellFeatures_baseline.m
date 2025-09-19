%% Load .mat files 

clear all;

% Set the directory for the experiment you need
fullCraniotomyBaselineDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\1._analysis.mat';
%fullCraniotomyCSDDir = 'V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\CSD\corrected_for_pinprick\1._analysis.mat';
% thinBoneBaselineDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\thinBone\_analysis.mat';

cd (fullCraniotomyBaselineDir);
% Get all .mat files in the directory
FilesAll = dir(fullfile(fullCraniotomyBaselineDir, '*_analysis.mat')); 

% Extract file names
fileNames = {FilesAll.name};
sortedFileNames = sortFileNames(fileNames);

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

%% after the extraction was done
clear all;
%for fullCraniotomy baseline data - load spreadsheet and add Var3
load('V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\AQuA2_data_fullCraniotomy_features_baseline.mat')
load('V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\multinucleated cells\multinucleatedCells.mat');
fullCraniotomy_combinedTable = addvars(combinedTable, multinucleatedCells.Var3, 'NewVariableNames', 'Multinucleated');
fullCraniotomy_multinucleated = fullCraniotomy_combinedTable{:,17};
fullCraniotomy_combinedTable_NM = fullCraniotomy_combinedTable(fullCraniotomy_multinucleated == 0, :);
fullCraniotomy_cellLocation = fullCraniotomy_combinedTable_NM{:,13};
fullCraniotomy_redLabel = fullCraniotomy_combinedTable_NM{:,14};
% By cell type
fullCraniotomy_perivascular_indices = fullCraniotomy_cellLocation == 0; % Indices of cells belonging to group 0
fullCraniotomy_nonPerivascular_indices = fullCraniotomy_cellLocation == 2; % Indices of cells belonging to group 2
fullCraniotomy_combinedTable_perivascular = fullCraniotomy_combinedTable(perivascular_indices,:);
fullCraniotomy_combinedTable_nonPerivascular = fullCraniotomy_combinedTable(nonPerivascular_indices,:);

%for thinBone data - create a new table 'multinucleatedCells' and add Var1
load('V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\thinBone\AQuA2_data_thinBone.mat')
thinBone_multinucleated = zeros(94,1);
thinBone_combinedTable = addvars(combinedTable, thinBone_multinucleated.Var1, 'NewVariableNames', 'Multinucleated');
thinBone_cellLocation = thinBone_combinedTable{:,13};
thinBone_all_NM_indices = thinBone_multinucleated == 0;
thinBone_numOnes = sum(thinBone_all_NM_indices);
thinBone_combinedTable = combinedTable(thinBone_all_NM_indices,:);

%% Count perivascular cells labeled with red and those not labeled with red
perivascular_red = sum(fullCraniotomy_cellLocation == 0 & fullCraniotomy_redLabel == 1);
perivascular_not_red = sum(fullCraniotomy_cellLocation == 0 & fullCraniotomy_redLabel == 0);
non_perivascular_red = sum(fullCraniotomy_cellLocation == 2 & fullCraniotomy_redLabel == 1);
non_perivascular_not_red = sum(fullCraniotomy_cellLocation == 2 & fullCraniotomy_redLabel == 0);

%% Features by Cell

%Perivascular
area_perivascular_NM = fullCraniotomy_combinedTable_perivascular(:,"Area(um2)");
perimeter_perivascular_NM = fullCraniotomy_combinedTable_perivascular(:,"Perimeter");
circularity_perivascular_NM = fullCraniotomy_combinedTable_perivascular(:,"Circularity");
maxDFF_perivascular_NM = fullCraniotomy_combinedTable_perivascular(:,"Max dFF");
dFFAUC_perivascular_NM = fullCraniotomy_combinedTable_perivascular(:,"dFF AUC");
duration_perivascular_NM = fullCraniotomy_combinedTable_perivascular(:,"Duration 10% to 10%");

%number of events 
numberOfEvents_perivascular_NM = fullCraniotomy_combinedTable_perivascular(:,"Number of Events");
totalEvents_perivascular_NM = sum(numberOfEvents_perivascular_NM{:,:});
%number of events in Hz
numberOfEvents_perivascular_Hz = table2cell(numberOfEvents_perivascular_NM);
numberOfEvents_perivascular_Hz = cell2mat(numberOfEvents_perivascular_Hz);
numberOfEvents_perivascular_Hz = numberOfEvents_perivascular_Hz / 900;
numberOfEvents_perivascular_Hz_new = numberOfEvents_perivascular_Hz * 1000;

%Non-Perivascular
area_nonPerivascular_NM = fullCraniotomy_combinedTable_nonPerivascular(:,"Area(um2)");
perimeter_nonPerivascular_NM = fullCraniotomy_combinedTable_nonPerivascular(:,"Perimeter");
circularity_nonPerivascular_NM = fullCraniotomy_combinedTable_nonPerivascular(:,"Circularity");
maxDFF_nonPerivascular_NM = fullCraniotomy_combinedTable_nonPerivascular(:,"Max dFF");
dFFAUC_nonPerivascular_NM = fullCraniotomy_combinedTable_nonPerivascular(:,"dFF AUC");
duration_nonPerivascular_NM = fullCraniotomy_combinedTable_nonPerivascular(:,"Duration 10% to 10%");

%number of events 
numberOfEvents_nonPerivascular_NM = fullCraniotomy_combinedTable_nonPerivascular(:,"Number of Events");
totalEvents_nonPerivascular_NM = sum(numberOfEvents_nonPerivascular_NM{:,:});
%number of events in Hz
numberOfEvents_nonPerivascular_Hz = table2cell(numberOfEvents_nonPerivascular_NM);
numberOfEvents_nonPerivascular_Hz = cell2mat(numberOfEvents_nonPerivascular_Hz);
numberOfEvents_nonPerivascular_Hz = numberOfEvents_nonPerivascular_Hz / 900;
numberOfEvents_nonPerivascular_Hz_new = numberOfEvents_nonPerivascular_Hz * 1000;


%% Perivascular By FOV

% Extract unique identifiers (e.g., up to "_run1")
perivascular_extractedNames = cellfun(@(x) extractBefore(x, '_r'), fullCraniotomy_combinedTable_perivascular.fileNameColumn, 'UniformOutput', false);
fullCraniotomy_combinedTable_perivascular.Group = perivascular_extractedNames;

perivascular_uniqueGroups = unique(fullCraniotomy_combinedTable_perivascular.Group);
perivascular_averagedTable = table();

for x = 1:length(perivascular_uniqueGroups)
    perivascular_groupName = perivascular_uniqueGroups{x};
    
    % Get rows corresponding to the current group
    perivascular_groupData = fullCraniotomy_combinedTable_perivascular(strcmp(fullCraniotomy_combinedTable_perivascular.Group, perivascular_groupName), :);
    perivascular_selectedVariableNames = perivascular_groupData.Properties.VariableNames(2:12);

    % Compute averages for numeric columns
     if size(perivascular_groupData, 1) == 1
        % Only one row, directly use the values
        perivascular_avgValues = perivascular_groupData{1, perivascular_selectedVariableNames};
     else
        % Compute averages for numeric columns
        perivascular_avgValues = mean(perivascular_groupData{:, perivascular_selectedVariableNames}, 1);
     end

    % Append to result table
    perivascular_tempTable = table({perivascular_groupName}, perivascular_avgValues(1),...
    perivascular_avgValues(2),perivascular_avgValues(3),perivascular_avgValues(4),...
    perivascular_avgValues(5),perivascular_avgValues(6),perivascular_avgValues(7),...
    perivascular_avgValues(8),perivascular_avgValues(9),perivascular_avgValues(10),...
    perivascular_avgValues(11),'VariableNames', [{'Group'}, perivascular_selectedVariableNames]);

    perivascular_averagedTable = [perivascular_averagedTable; perivascular_tempTable];
end

%% nonPerivascular By FOV

% Extract unique identifiers (e.g., up to "_run1")
nonPerivascular_extractedNames = cellfun(@(x) extractBefore(x, '_r'), fullCraniotomy_combinedTable_nonPerivascular.fileNameColumn, 'UniformOutput', false);
fullCraniotomy_combinedTable_nonPerivascular.Group = nonPerivascular_extractedNames;

nonPerivascular_uniqueGroups = unique(fullCraniotomy_combinedTable_nonPerivascular.Group);
nonPerivascular_averagedTable = table();

for x = 1:length(nonPerivascular_uniqueGroups)
    nonPerivascular_groupName = nonPerivascular_uniqueGroups{x};
    
    % Get rows corresponding to the current group
    nonPerivascular_groupData = fullCraniotomy_combinedTable_nonPerivascular(strcmp(fullCraniotomy_combinedTable_nonPerivascular.Group, nonPerivascular_groupName), :);
    nonPerivascular_selectedVariableNames = nonPerivascular_groupData.Properties.VariableNames(2:12);

    % Compute averages for numeric columns
     if size(nonPerivascular_groupData, 1) == 1
        % Only one row, directly use the values
        nonPerivascular_avgValues = nonPerivascular_groupData{1, nonPerivascular_selectedVariableNames};
     else
        % Compute averages for numeric columns
        nonPerivascular_avgValues = mean(nonPerivascular_groupData{:, nonPerivascular_selectedVariableNames}, 1);
     end

    % Append to result table
    nonPerivascular_tempTable = table({nonPerivascular_groupName}, nonPerivascular_avgValues(1), ...
        nonPerivascular_avgValues(2),nonPerivascular_avgValues(3),nonPerivascular_avgValues(4),...
    nonPerivascular_avgValues(5),nonPerivascular_avgValues(6),nonPerivascular_avgValues(7),...
    nonPerivascular_avgValues(8),nonPerivascular_avgValues(9),nonPerivascular_avgValues(10),...
    nonPerivascular_avgValues(11),'VariableNames', [{'Group'}, nonPerivascular_selectedVariableNames]);

    nonPerivascular_averagedTable = [nonPerivascular_averagedTable; nonPerivascular_tempTable];
end

%% Features by FOV

% Perivascular
area_perivascular = perivascular_averagedTable{:,2};
perimeter_perivascular = perivascular_averagedTable{:,3};
circularity_perivascular = perivascular_averagedTable{:,4};
maxDFF_perivascular = perivascular_averagedTable{:,5};
dFFAUC_perivascular = perivascular_averagedTable{:,11};
numberOfEvents_perivascular = perivascular_averagedTable{:,12};

% nonPerivascular
area_nonPerivascular = nonPerivascular_averagedTable{:,2};
perimeter_nonPerivascular = nonPerivascular_averagedTable{:,3};
circularity_nonPerivascular = nonPerivascular_averagedTable{:,4};
maxDFF_nonPerivascular = nonPerivascular_averagedTable{:,5};
dFFAUC_nonPerivascular = nonPerivascular_averagedTable{:,11};
numberOfEvents_nonPerivascular = nonPerivascular_averagedTable{:,12};


%% Features of perivascular cells (positive vs negative coefficient)

% Read the Excel file into a table
coefficients = readtable('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\coefficients\coefficients.xlsx');
save('coefficients.mat', 'data');

load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\coefficients\coefficients.mat');
fullCraniotomy_combinedTable2 = addvars(fullCraniotomy_combinedTable, coefficients.Var3, 'NewVariableNames', 'Coefficients(1=pos;2=neg)'); % 1 = positive, 2 = negative

fullCraniotomy_coefficients = fullCraniotomy_combinedTable2{:,18};
perivascular_positive_indices = fullCraniotomy_coefficients == 1;
perivascular_negative_indices = fullCraniotomy_coefficients == 2;

fullCraniotomy_combinedTable2_perivascular_positive_indices = fullCraniotomy_combinedTable2(perivascular_positive_indices,:);
fullCraniotomy_combinedTable2_perivascular_negative_indices = fullCraniotomy_combinedTable2(perivascular_negative_indices,:);


%Perivascular Positive
area_perivascular_positive= fullCraniotomy_combinedTable2_perivascular_positive_indices(:,"Area(um2)");
perimeter_perivascular_positive = fullCraniotomy_combinedTable2_perivascular_positive_indices(:,"Perimeter");
circularity_perivascular_positive = fullCraniotomy_combinedTable2_perivascular_positive_indices(:,"Circularity");
numberOfEvents_perivascular_positive = fullCraniotomy_combinedTable2_perivascular_positive_indices(:,"Number of Events");
maxDFF_perivascular_positive = fullCraniotomy_combinedTable2_perivascular_positive_indices(:,"Max dFF");
dFFAUC_perivascular_positive = fullCraniotomy_combinedTable2_perivascular_positive_indices(:,"dFF AUC");

%Perivascular Negative
area_perivascular_negative = fullCraniotomy_combinedTable2_perivascular_negative_indices(:,"Area(um2)");
perimeter_perivascular_negative = fullCraniotomy_combinedTable2_perivascular_negative_indices(:,"Perimeter");
circularity_perivascular_negative = fullCraniotomy_combinedTable2_perivascular_negative_indices(:,"Circularity");
numberOfEvents_perivascular_negative = fullCraniotomy_combinedTable2_perivascular_negative_indices(:,"Number of Events");
maxDFF_perivascular_negative = fullCraniotomy_combinedTable2_perivascular_negative_indices(:,"Max dFF");
dFFAUC_perivascular_negative = fullCraniotomy_combinedTable2_perivascular_negative_indices(:,"dFF AUC");

%number of events in Hz
%positive
numberOfEvents_perivascular_positive_Hz = table2cell(numberOfEvents_perivascular_positive);
numberOfEvents_perivascular_positive_Hz = cell2mat(numberOfEvents_perivascular_positive_Hz);
numberOfEvents_perivascular_positive_Hz = numberOfEvents_perivascular_positive_Hz / 900;
%negative
numberOfEvents_perivascular_negative_Hz = table2cell(numberOfEvents_perivascular_negative);
numberOfEvents_perivascular_negative_Hz = cell2mat(numberOfEvents_perivascular_negative_Hz);
numberOfEvents_perivascular_negative_Hz = numberOfEvents_perivascular_negative_Hz / 900;


%% duplicating table to manually delete a cell if necessary

figure;
plot(combinedTable.dFF{425,:})
combinedTable2 = combinedTable; %duplicating table to manually delete a cell if necessary

%% Sort combinedTable by number of events

combinedTable_sorted = sortrows(combinedTable, 'Cell location (0,perivascular;1,adjacent;2,none)', 'descend');

%% Filter data by cell location

varTypes = varfun(@class, combinedTable_sorted, 'OutputFormat', 'cell'); %check data types

%create table with all variables but cell ID and dFF
dataAll = combinedTable_sorted{:, 2:14};
dataAll_Zscore = zscore(dataAll);
dataAll2 = combinedTable_sorted{:, [2:9,11:14]};


% Separate cells by location
cellLocation = combinedTable_sorted{:, 13};

logicalIndex0 = cellLocation == 0; %perivascular
rowsWithValue0 = combinedTable_sorted(logicalIndex0, :);
dataPerivascular = rowsWithValue0{:, 2:14};

logicalIndex2 = cellLocation == 2; %non-perivascular
rowsWithValue2 = combinedTable_sorted(logicalIndex2, :);
dataNonPerivascular = rowsWithValue2{:, 2:14};

%% Renumber cells in combinedTable to simplify cell identification

%create cellID to substitute cell numbers
cellID_all = (1:544)';
cellID_all = num2cell(cellID_all);
cellID = combinedTable_sorted{:, 1};
cellID_combined = [cellID_all, cellID];
cellID_combinedTable = cell2table(cellID_combined);
