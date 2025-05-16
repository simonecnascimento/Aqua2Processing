% Input data for each category: Chronic [up, no change, down]
P_acute_up = [6, 9, 6];
P_acute_no = [7, 0, 2];
P_acute_down = [3, 9, 22];

% Combine data into a matrix
P_matrix = [P_acute_up; P_acute_no; P_acute_down];

% Calculate the total number for each acute condition
P_total_up = sum(P_acute_up);
P_total_no = sum(P_acute_no);
P_total_down = sum(P_acute_down);

% Conditional probability calculation
P_conditional = (P_matrix ./ sum(P_matrix, 2))*100;

% Input data for each category: Chronic [up, no change, down]
NP_acute_up = [2, 26, 4];
NP_acute_no = [26, 1, 5];
NP_acute_down = [6, 12, 103];

% Combine data into a matrix
NP_matrix = [NP_acute_up; NP_acute_no; NP_acute_down];

% Calculate the total number for each acute condition
NP_total_up = sum(NP_acute_up);
NP_total_no = sum(NP_acute_no);
NP_total_down = sum(NP_acute_down);

% Conditional probability calculation
NP_conditional = (NP_matrix ./ sum(NP_matrix, 2))*100;

% Define row and column labels
row_labels = {'↑', 'x', '↓'};  % Acute conditions
col_labels = {'↑', 'x', '↓'};  % Chronic conditions

% Find the global min and max values across both matrices for consistent color scaling
min_val = min([P_conditional(:); NP_conditional(:)]);
max_val = max([P_conditional(:); NP_conditional(:)]);

% Plot heatmaps for visual comparison
figure;

subplot(1,2,1);
h1 = heatmap(col_labels, row_labels, P_conditional, 'Title', 'Perivascular', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
h1.ColorLimits = [min_val, max_val];  % Set the same color limits for both matrices

subplot(1,2,2);
h2 = heatmap(col_labels, row_labels, NP_conditional, 'Title', 'Non-perivascular', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
h2.ColorLimits = [min_val, max_val];  % Set the same color limits for both matrices


%%

% Input data for each category: Chronic [up, no change, down]
All_acute_up = [8, 35, 10];
All_acute_no = [33, 1, 7];
All_acute_down = [9, 21, 125];

% Combine data into a matrix
All_matrix = [All_acute_up; All_acute_no; All_acute_down];

% Calculate the total number for each acute condition
All_total_up = sum(All_acute_up);
All_total_no = sum(All_acute_no);
All_total_down = sum(All_acute_down);

% Conditional probability calculation
All_conditional = (All_matrix ./ sum(All_matrix, 2))*100;

% Define row and column labels
row_labels = {'↑', 'x', '↓'};  % Acute conditions
col_labels = {'↑', 'x', '↓'};  % Chronic conditions

% Find the global min and max values across both matrices for consistent color scaling
min_val = min(All_conditional(:));
max_val = max(All_conditional(:));

% Plot heatmaps for visual comparison
figure;

h1 = heatmap(col_labels, row_labels, All_conditional, 'Title', 'All cells', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
h1.ColorLimits = [min_val, max_val];  % Set the same color limits for both matrices


