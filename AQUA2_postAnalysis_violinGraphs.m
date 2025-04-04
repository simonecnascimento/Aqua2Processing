%% Data analysis fullCraniotomy vs thinBone

clear all; 
%for fullCraniotomy data - load spreadsheet and add Var3
load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\AQuA2_data_fullCraniotomy_features_baseline.mat')
load('D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\multinucleated cells\multinucleatedCells.mat');
fullCraniotomy_combinedTable = addvars(combinedTable, multinucleatedCells.Var3, 'NewVariableNames', 'Multinucleated');
fullCraniotomy_multinucleated = fullCraniotomy_combinedTable{:,17};
fullCraniotomy_cellLocation = fullCraniotomy_combinedTable{:,13};
fullCraniotomy_all_NM_indices = fullCraniotomy_multinucleated == 0;
fullCraniotomy_numOnes = sum(fullCraniotomy_all_NM_indices);


%% Circularity Perivascular vs Non-perivascular

circularity = combinedTable{:,4};
circularity_perivascular_NM = circularity(perivascular_NM_indices); % Data points belonging to group 0
circularity_nonPerivascular_NM = circularity(nonPerivascular_NM_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin({circularity_perivascular_NM, circularity_nonPerivascular_NM}, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
ylabel('Circularity', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 0.7]);
xticks([1, 2]);
xticklabels({'Perivascular', 'Non-Perivascular'});

% Remove top and right lines
box off;

% Perform t-test
[h, circularity_p_value] = ttest2(circularity_perivascular_NM, circularity_nonPerivascular_NM);

%% Circularity Non-perivascular ONLY

circularity = combinedTable_fullCraniotomy{:,4};
circularity_perivascular = circularity(perivascular_indices); % Data points belonging to group 0
circularity_nonPerivascular = circularity(nonPerivascular_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin(circularity_nonPerivascular, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
%ylabel('Circularity', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 0.7]);
xticks([1, 2]);
xticklabels('');
title('Circularity', FontSize=30);

% Remove top and right lines
box off;

% Perform t-test
[h, circularity_p_value] = ttest2(circularity_perivascular, circularity_nonPerivascular);

%% Area

area = combinedTable_fullCraniotomy{:,2};
area_perivascular_NM = area(perivascular_NM_indices); % Data points belonging to group 0
area_nonPerivascular_NM = area(nonPerivascular_NM_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin({area_perivascular_NM, area_nonPerivascular_NM}, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
ylabel('Area (um2)', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 1000]);
xticks([1, 2]);
xticklabels({'Perivascular', 'Non-Perivascular'});

% Remove top and right lines
box off;

% Perform t-test
[h, area_p_value] = ttest2(area_perivascular_NM, area_nonPerivascular_NM);

%% Area non perivascular ONLY

area = combinedTable_fullCraniotomy{:,2};
area_perivascular = area(perivascular_indices); % Data points belonging to group 0
area_nonPerivascular = area(nonPerivascular_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin(area_nonPerivascular, 'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
%ylabel('Area (um2)', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 1000]);
xticks([1, 2]);
%xticklabels({'Perivascular', 'Non-Perivascular'});
title('Area(um2)', FontSize=30);

% Remove top and right lines
box off;

% Adjust y-axis tick font size and remove x-axis ticks
ax = gca;
ax.YAxis.FontSize = 30;  % Set y-axis number font size to 30
set(gca, 'xtick', []);   % Remove x-axis ticks

% Perform t-test
[h, area_p_value] = ttest2(area_perivascular, area_nonPerivascular);


%% Perimeter

perimeter = combinedTable_fullCraniotomy{:,3};
perimeter_perivascular_NM = perimeter(perivascular_NM_indices); % Data points belonging to group 0
perimeter_nonPerivascular_NM = perimeter(nonPerivascular_NM_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin({perimeter_perivascular_NM, perimeter_nonPerivascular_NM}, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
ylabel('Perimeter', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 500]);
xticks([1, 2]);
xticklabels({'Perivascular', 'Non-Perivascular'});

% Remove top and right lines
box off;

% Perform t-test
[h, perimeter_p_value] = ttest2(perimeter_perivascular_NM, perimeter_nonPerivascular_NM);

%% Perimeter non perivascular ONLY

perimeter = combinedTable_fullCraniotomy{:,3};
perimeter_perivascular = perimeter(perivascular_indices); % Data points belonging to group 0
perimeter_nonPerivascular = perimeter(nonPerivascular_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin(perimeter_nonPerivascular, 'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
%ylabel('Perimeter', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 500]);
xticks([1, 2]);
%xticklabels({'Perivascular', 'Non-Perivascular'});
title('Perimeter(um)', FontSize=30);

% Remove top and right lines
box off;

% Adjust y-axis tick font size and remove x-axis ticks
ax = gca;
ax.YAxis.FontSize = 30;  % Set y-axis number font size to 30
set(gca, 'xtick', []);   % Remove x-axis ticks

% Perform t-test
[h, perimeter_p_value] = ttest2(perimeter_perivascular, perimeter_nonPerivascular);

%% Number of Events

numberOfEvents = combinedTable_fullCraniotomy{:,12};
numberOfEvents_perivascular_NM = numberOfEvents(perivascular_NM_indices); % Data points belonging to group 0
numberOfEvents_nonPerivascular_NM = numberOfEvents(nonPerivascular_NM_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin({numberOfEvents_perivascular_NM, numberOfEvents_nonPerivascular_NM}, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
ylabel('Number of Events', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 50]);
xticks([1, 2]);
xticklabels({'Perivascular', 'Non-Perivascular'});

% Remove top and right lines
box off;

% Perform t-test
[h, numberOfEvents_p_value] = ttest2(numberOfEvents_perivascular_NM, numberOfEvents_nonPerivascular_NM);

%% Number of Events non perivascular ONLY

numberOfEvents = combinedTable_fullCraniotomy{:,12};
numberOfEvents_perivascular = numberOfEvents(perivascular_indices); % Data points belonging to group 0
numberOfEvents_nonPerivascular = numberOfEvents(nonPerivascular_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin(numberOfEvents_nonPerivascular, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
%ylabel('Number of Events', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 30]);
xticks([1, 2]);
%xticklabels({'Perivascular', 'Non-Perivascular'});
%xticklabels({'Perivascular', 'Non-Perivascular'});
title('Number of events/15min', FontSize=30);

% Remove top and right lines
box off;

% Adjust y-axis tick font size and remove x-axis ticks
ax = gca;
ax.YAxis.FontSize = 30;  % Set y-axis number font size to 30
set(gca, 'xtick', []);   % Remove x-axis ticks

% Remove top and right lines
box off;

% Perform t-test
[h, numberOfEvents_p_value] = ttest2(numberOfEvents_perivascular, numberOfEvents_nonPerivascular);

%% Number of Events non perivascular ONLY (in Hz)

numberOfEvents = combinedTable_fullCraniotomy{:,12};
numberOfEvents_perivascular = numberOfEvents(perivascular_indices); % Data points belonging to group 0
%numberOfEvents_nonPerivascular = numberOfEvents(nonPerivascular_indices); % Data points belonging to group 2

% Convert number of events per 15 minutes to frequency (Hz)
numberOfEvents_perivascular_Hz = numberOfEvents_perivascular / 900; % Convert to Hz
numberOfEvents_nonPerivascular_Hz = numberOfEvents_nonPerivascular / 900; % Convert to Hz


% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin(numberOfEvents_nonPerivascular_Hz, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
%ylabel('Number of Events', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 30/900]);
xticks([1, 2]);
%xticklabels({'Perivascular', 'Non-Perivascular'});
%xticklabels({'Perivascular', 'Non-Perivascular'});
title('Number of events/15min', FontSize=30);

% Remove top and right lines
box off;

% Adjust y-axis tick font size and remove x-axis ticks
ax = gca;
ax.YAxis.FontSize = 30;  % Set y-axis number font size to 30
set(gca, 'xtick', []);   % Remove x-axis ticks

% Remove top and right lines
box off;

% Perform t-test
[h, numberOfEvents_p_value] = ttest2(numberOfEvents_perivascular, numberOfEvents_nonPerivascular);

%% Max DFF

maxDFF = combinedTable{:,5};
maxDFF_perivascular_NM = maxDFF(perivascular_NM_indices); % Data points belonging to group 0
maxDFF_nonPerivascular_NM = maxDFF(nonPerivascular_NM_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin({maxDFF_perivascular_NM, maxDFF_nonPerivascular_NM}, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
ylabel('Max dFF', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 5]);
xticks([1, 2]);
xticklabels({'Perivascular', 'Non-Perivascular'});

% Remove top and right lines
box off;

% Perform t-test
[h, maxDFF_p_value] = ttest2(maxDFF_perivascular_NM, maxDFF_nonPerivascular_NM);

%% Max DFF non perivascular ONLY

maxDFF = combinedTable{:,5};
maxDFF_perivascular = maxDFF(perivascular_indices); % Data points belonging to group 0
maxDFF_nonPerivascular = maxDFF(nonPerivascular_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin(maxDFF_nonPerivascular, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
%ylabel('Max dFF', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 5]);
xticks([1, 2]);
%xticklabels({'Perivascular', 'Non-Perivascular'});
title('Max fluorescence', FontSize=30); %title('Max dFF', FontSize=30);

% Adjust y-axis tick font size and remove x-axis ticks
ax = gca;
ax.YAxis.FontSize = 30;  % Set y-axis number font size to 30
set(gca, 'xtick', []);   % Remove x-axis ticks

% Remove top and right lines
box off;


% Perform t-test
[h, maxDFF_p_value] = ttest2(maxDFF_perivascular, maxDFF_nonPerivascular);


%% dFF AUC

dFF_AUC = combinedTable{:,11};
dFF_AUC_perivascular_NM = dFF_AUC(perivascular_NM_indices); % Data points belonging to group 0
dFF_AUC_nonPerivascular_NM = dFF_AUC(nonPerivascular_NM_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin({dFF_AUC_perivascular_NM, dFF_AUC_nonPerivascular_NM}, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
ylabel('dFF AUC', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 250]);
xticks([1, 2]);
xticklabels({'Perivascular', 'Non-Perivascular'});

% Remove top and right lines
box off;

% Perform t-test
[h, dFF_AUC_p_value] = ttest2(dFF_AUC_perivascular_NM, dFF_AUC_nonPerivascular_NM);

%% dFF AUC non perivascular ONLY

dFF_AUC = combinedTable{:,11};
dFF_AUC_perivascular = dFF_AUC(perivascular_indices); % Data points belonging to group 0
dFF_AUC_nonPerivascular = dFF_AUC(nonPerivascular_indices); % Data points belonging to group 2

% Create a combined violin plot for both groups
figure;
barColors = [1 0 0; 0 0 1];
violin({dFF_AUC_perivascular, dFF_AUC_nonPerivascular}, ...
    'facecolor', barColors, 'edgecolor', [], 'mc', 'k');
ylabel('dFF AUC', 'FontName', 'Helvetica', 'FontSize', 18);
ylim([0, 250]);
xticks([1, 2]);
xticklabels({'Perivascular', 'Non-Perivascular'});

title('Max dFF', FontSize=30);

% Adjust y-axis tick font size and remove x-axis ticks
ax = gca;
ax.YAxis.FontSize = 30;  % Set y-axis number font size to 30
set(gca, 'xtick', []);   % Remove x-axis ticks

% Remove top and right lines
box off;

% Perform t-test
[h, dFF_AUC_p_value] = ttest2(dFF_AUC_perivascular, dFF_AUC_nonPerivascular);