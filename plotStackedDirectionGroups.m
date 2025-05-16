function plotStackedDirectionGroups(directionCounts, savePath)
    % plotStackedDirectionGroups - Plot stacked bar chart of direction counts grouped in 3s
    %
    % Inputs:
    %   directionCounts - a matrix of size (rows x 9 columns), where rows are samples and
    %                     each column corresponds to a direction pair (e.g., ↑↑, ↑x, ↑↓, etc.)

    % Sum the rows to get combined counts
    totalData = sum(directionCounts, 1);  % 1 row with 9 columns

    % Define 3 groups of 3 columns
    groups = {[1 2 3], [4 5 6], [7 8 9]};
    groupLabels = {'↑', 'x', '↓'};  % Group by first arrow
    subgroupLabels = {'↑', 'x', '↓'};  % Second arrow for legend

    % Compute percentage within each group
    groupedPercents = zeros(3, 3);  % 3 groups (columns), 3 elements each (rows)

    for g = 1:3
        cols = groups{g};
        values = totalData(cols)';
        total = sum(values);
        if total == 0
            groupedPercents(:, g) = 0;
        else
            groupedPercents(:, g) = ((values / total) * 100);
        end
    end

    % Plot stacked bar chart
    fig = figure;
    bar(groupedPercents', 'stacked');
    ylabel('Percentage');
    ylim([0,120]);
    xticklabels(groupLabels);
    legend(subgroupLabels, 'Location', 'eastoutside');
    title('Stacked Bar: Direction Group Percentages');

    if nargin >= 2 && ~isempty(savePath)
        if ~exist(savePath, 'dir')
            mkdir(savePath);
        end
        filename = fullfile(savePath, 'Stacked Bar - P and NP together.eps');
        saveas(fig, filename);
        fprintf('Figure saved to: %s\n', filename);
    end
end
