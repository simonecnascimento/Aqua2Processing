function [All_matrix, P_matrix, NP_matrix] = plotCluster6Distribution(clusterCounts, outputDir)

    % Define row and column labels (acute x chronic)
    % directions = ["↑↑", "↑x", "↑↓", "x↑", "xx", "x↓"];
    row_labels = {'↑', 'X'};
    col_labels = {'↑', 'X', '↓'};

    % Split cluster counts
    % Order: ["↑↑", "↑x", "↑↓", "x↑", "xx", "x↓"];
    P_counts  = clusterCounts(1, :); % perivascular
    NP_counts = clusterCounts(2, :); % non-perivascular

    % Format into 2x3 matrices: rows = acute, cols = chronic
    P_matrix  = reshape(P_counts,  [3, 2])'; % reshape as (chronic cols × acute rows), then transpose
    NP_matrix = reshape(NP_counts, [3, 2])';

    % ------------------ CONDITIONAL PROBABILITIES ------------------
    P_conditional  = round((P_matrix  ./ sum(P_matrix,  2)) * 100, 1);
    NP_conditional = round((NP_matrix ./ sum(NP_matrix, 2)) * 100, 1);

    % ------------------ OVERALL PERCENTAGE ------------------
    P_percent  = round((P_matrix  / sum(P_matrix(:)))  * 100, 1);
    NP_percent = round((NP_matrix / sum(NP_matrix(:))) * 100, 1);

    % ------------------ GLOBAL COLOR SCALE ------------------
    min_cond = min([P_conditional(:); NP_conditional(:)]);
    max_cond = max([P_conditional(:); NP_conditional(:)]);
    
    min_perc = 0; %min([P_percent(:); NP_percent(:)]);
    max_perc = 70; %max([P_percent(:); NP_percent(:)]);

    % ------------------ HEATMAPS ------------------
    % Conditional - Perivascular
    figP = figure;
    hP = heatmap(col_labels, row_labels, P_conditional, ...
        'Title', 'Perivascular (conditional %)', ...
        'XLabel', 'Chronic', 'YLabel', 'Acute');
    hP.ColorLimits = [min_cond, max_cond];
    saveas(figP, fullfile(outputDir, 'Heatmap_ConditionalP_6clusters.eps'), 'epsc');

    % Conditional - Non-perivascular
    figNP = figure;
    hNP = heatmap(col_labels, row_labels, NP_conditional, ...
        'Title', 'Non-perivascular (conditional %)', ...
        'XLabel', 'Chronic', 'YLabel', 'Acute');
    hNP.ColorLimits = [min_cond, max_cond];
    saveas(figNP, fullfile(outputDir, 'Heatmap_ConditionalNP_6clusters.eps'), 'epsc');

    % Percentage - Perivascular
    figP2 = figure;
    hP2 = heatmap(col_labels, row_labels, P_percent, ...
        'Title', 'Perivascular (total %)', ...
        'XLabel', 'Chronic', 'YLabel', 'Acute');
    hP2.ColorLimits = [min_perc, max_perc];
    saveas(figP2, fullfile(outputDir, 'Heatmap_PercentageP_6clusters.eps'), 'epsc');

    % Percentage - Non-perivascular
    figNP2 = figure;
    hNP2 = heatmap(col_labels, row_labels, NP_percent, ...
        'Title', 'Non-perivascular (total %)', ...
        'XLabel', 'Chronic', 'YLabel', 'Acute');
    hNP2.ColorLimits = [min_perc, max_perc];
    saveas(figNP2, fullfile(outputDir, 'Heatmap_PercentageNP_6clusters.eps'), 'epsc');

    % ------------------ COMBINED ------------------
    All_matrix = P_matrix + NP_matrix;
    All_conditional = round((All_matrix ./ sum(All_matrix, 2)) * 100, 1);
    All_percent = round((All_matrix / sum(All_matrix(:))) * 100, 1);

    % Heatmaps for All
    figAllC = figure;
    hAllC = heatmap(col_labels, row_labels, All_conditional, ...
        'Title', 'All Cells (conditional %)', ...
        'XLabel', 'Chronic', 'YLabel', 'Acute');
    hAllC.ColorLimits = [min_cond, max_cond];
    saveas(figAllC, fullfile(outputDir, 'Heatmap_ConditionalAll_6clusters.eps'), 'epsc');

    figAllP = figure;
    hAllP = heatmap(col_labels, row_labels, All_percent, ...
        'Title', 'All Cells (total %)', ...
        'XLabel', 'Chronic', 'YLabel', 'Acute');
    hAllP.ColorLimits = [min_perc, max_perc];
    saveas(figAllP, fullfile(outputDir, 'Heatmap_PercentageAll_6clusters.eps'), 'epsc');

    % ------------------ Residuals ------------------
    % Optional residual plot — uncomment if needed
    % plotStdResiduals(NP_matrix, row_labels, col_labels, outputDir);

end
