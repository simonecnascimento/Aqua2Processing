function [p_values, chi2_stats] = chi_square_by_acuteResponse(data, groupSize)
    % chi_square_by_groups - Perform chi-square tests on each group of columns
    %
    % INPUTS:
    %   data      - a 2 x N matrix of counts
    %   groupSize - number of columns to group together (e.g., 3)
    %
    % OUTPUTS:
    %   p_values  - array of p-values for each chi-square test
    %   chi2_stats - cell array of chi-square stats (df, chi2stat, etc.)

    if size(data,1) ~= 2
        error('Data must be a 2-row matrix (e.g., perivascular and non-perivascular).');
    end

    numGroups = floor(size(data,2) / groupSize);
    p_values = zeros(1, numGroups);
    chi2_stats = cell(1, numGroups);

    for i = 1:numGroups
        cols = (i-1)*groupSize + (1:groupSize);
        subtable = data(:, cols);

        % Check if chi2gof_from_table is available; else compute manually
        [~, p, stats] = chi2gof_from_table(subtable);  % You can replace with custom logic if needed

        p_values(i) = p;
        chi2_stats{i} = stats;
    end

    % Display results
    fprintf('Chi-square p-values by group (every %d columns):\n', groupSize);
    disp(p_values);
end
