function [p_values, chi2_stats] = chi_square_by_chronicResponse(data)
    % chi_square_predefined_groups - Perform chi-square tests on predefined 3-column groups
    %
    % INPUT:
    %   data - a 2 x 9 matrix of counts
    %
    % OUTPUT:
    %   p_values   - array of p-values for each chi-square test
    %   chi2_stats - cell array of chi-square stats for each group

    if size(data, 1) ~= 2 || size(data, 2) ~= 9
        error('Input must be a 2x9 matrix of counts.');
    end

    % Define internal grouping: group 1 = cols 1,4,7; group 2 = 2,5,8; group 3 = 3,6,9
    groupCols = {
        [1 4 7];
        [2 5 8];
        [3 6 9]
    };

    numGroups = numel(groupCols);
    p_values = zeros(1, numGroups);
    chi2_stats = cell(1, numGroups);

    for i = 1:numGroups
        cols = groupCols{i};
        subtable = data(:, cols);

        % Perform chi-square test of independence
        [~, p, stats] = chi2gof_from_table(subtable);

        p_values(i) = p;
        chi2_stats{i} = stats;
    end

    % Print results
    fprintf('Chi-square p-values for predefined groups:\n');
    disp(p_values);
end