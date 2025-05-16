function [All_matrix, P_matrix, NP_matrix] = plotConditionalProbabilities(directionCounts, outputDir)
    
    % Define row and column labels
    row_labels = {'↑', 'x', '↓'};  % Acute conditions
    col_labels = {'↑', 'x', '↓'};  % Chronic conditions

    % Perivascular
    % Input data for each category: Acute (up, no change, down)
    P_acute_up = directionCounts(1, 1:3); % First row: perivascular up
    P_acute_no = directionCounts(1, 4:6); % First row: perivascular no change
    P_acute_down = directionCounts(1, 7:9); % First row: perivascular down

    % Combine data into a matrix for perivascular
    P_matrix = [P_acute_up; P_acute_no; P_acute_down];
    
    % Conditional probability calculation for perivascular cells
    P_conditional = round((P_matrix ./ sum(P_matrix, 2)) * 100, 1);
    P_percent = round((P_matrix / sum(P_matrix(:))) * 100, 1);

    % Non-perivascular
    % Input data for each category: Non-perivascular (up, no change, down)
    NP_acute_up = directionCounts(2, 1:3); % Second row: non-perivascular up
    NP_acute_no = directionCounts(2, 4:6); % Second row: non-perivascular no change
    NP_acute_down = directionCounts(2, 7:9); % Second row: non-perivascular down

    % Combine data into a matrix for non-perivascular
    NP_matrix = [NP_acute_up; NP_acute_no; NP_acute_down];

    % Conditional probability calculation for non-perivascular cells
    NP_conditional = round((NP_matrix ./ sum(NP_matrix, 2)) * 100, 1);
    NP_percent = round((NP_matrix / sum(NP_matrix(:))) * 100, 1);

            % Plot conditionals
            % Find the global min and max values across both matrices for consistent color scaling
            min_conditional = min([P_conditional(:); NP_conditional(:)]);
            max_conditional = max([P_conditional(:); NP_conditional(:)]);
            
            % Plot heatmaps for visual comparison
            figConditionalP = figure;
        
            % Plot perivascular heatmap
            %subplot(1, 2, 1);
            h1 = heatmap(col_labels, row_labels, P_conditional, 'Title', 'Perivascular', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
            h1.ColorLimits = [min_conditional, max_conditional];  % Set the same color limits for both matrices
            saveas(figConditionalP, 'Heatmap_ConditionalP', 'epsc');
        
            figConditionalNP = figure;
            % Plot non-perivascular heatmap
            %subplot(1, 2, 2);
            h2 = heatmap(col_labels, row_labels, NP_conditional, 'Title', 'Non-perivascular', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
            h2.ColorLimits = [min_conditional, max_conditional];  % Set the same color limits for both matrices
            saveas(figConditionalNP, 'Heatmap_ConditionalNP', 'epsc');

            % Plot percentages
            % Find the global min and max values across both matrices for consistent color scaling
            min_percentage = min([P_percent(:); NP_percent(:)]);
            max_percentage = max([P_percent(:); NP_percent(:)]);
            
            % Plot heatmaps for visual comparison
            figPercentagesP = figure;
        
            % Plot perivascular heatmap
            %subplot(1, 2, 1);
            h1 = heatmap(col_labels, row_labels, P_percent, 'Title', 'Perivascular', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
            h1.ColorLimits = [min_percentage, max_percentage];  % Set the same color limits for both matrices
            saveas(figPercentagesP, 'Heatmap_PercentageP', 'epsc');
        
            figPercentagesNP = figure;
            % Plot non-perivascular heatmap
            %subplot(1, 2, 2);
            h2 = heatmap(col_labels, row_labels, NP_percent, 'Title', 'Non-perivascular', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
            h2.ColorLimits = [min_percentage, max_percentage];  % Set the same color limits for both matrices
            saveas(figPercentagesNP, 'Heatmap_PercentageNP', 'epsc');

    % Input data for each category: Chronic [up, no change, down]
    % element‑wise sums (perivascular + non‑perivascular)
    All_acute_up   = P_acute_up   + NP_acute_up;    % [6 9 6] + [2 26 4]   -> [8 35 10]
    All_acute_no   = P_acute_no   + NP_acute_no;    % use the matching “no‑change” rows
    All_acute_down = P_acute_down + NP_acute_down;  % use the matching “down” rows
    
    % Combine data into a matrix
    All_matrix = [All_acute_up; All_acute_no; All_acute_down];
    
    % Calculate the total number for each acute condition
    All_total_up = sum(All_acute_up);
    All_total_no = sum(All_acute_no);
    All_total_down = sum(All_acute_down);
    
    % Conditional probability calculation
    All_conditional = round((All_matrix ./ sum(All_matrix, 2)) * 100, 1);
    All_percent = round((All_matrix / sum(All_matrix(:))) * 100, 1);

            % Plot conditionals heatmap
            min_conditionalAll = min(All_conditional(:));
            max_conditionalAll = max(All_conditional(:));
            figConditionalAll = figure;
            h1 = heatmap(col_labels, row_labels, All_conditional, 'Title', 'All cells', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
            h1.ColorLimits = [min_conditionalAll, max_conditionalAll];  % Set the same color limits for both matrices
            saveas(figConditionalAll, 'Heatmap_ConditionalAll', 'epsc');

            % Plot percentages heatmpa
            min_percentageAll = min(All_percent(:));
            max_percentageAll = max(All_percent(:));
            figPercentagesAll = figure;
            h1 = heatmap(col_labels, row_labels, All_percent, 'Title', 'All cells', 'XLabel', 'Chronic Response', 'YLabel', 'Acute Response');
            h1.ColorLimits = [min_percentage, max_percentage];  % Set the same color limits for both matrices
            saveas(figPercentagesAll, 'Heatmap_PercentageAll2', 'epsc');

    % Calculate and plot residuals
    plotStdResiduals(NP_matrix, row_labels, col_labels, outputDir);

end
