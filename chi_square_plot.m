function [p_values, chi2_stats] = chi_square_plot(directionCounts)
    % Initialize variables
    p_values = zeros(1, size(directionCounts, 2));
    chi2_stats = cell(1, size(directionCounts, 2));
    
    % Compute row and column sums for expected counts
    rowSums = sum(directionCounts, 2);  % Sum of rows
    colSums = sum(directionCounts, 1);  % Sum of columns
    total = sum(directionCounts(:));    % Grand total
    
    % Compute expected counts matrix (for each cell)
    expected = (rowSums * colSums) / total;
    
    % Perform Chi-Square Goodness of Fit for each cluster (column)
    for i = 1:size(directionCounts, 2)
        observed = directionCounts(:, i);  % get counts for cluster i
        
        % Only test if there are at least two non-zero counts
        if sum(observed) > 0 && sum(observed > 0) > 1
            % Perform chi-square test using expected values for each column
            [~, p, stats] = chi2gof([1, 2], 'freq', observed, 'expected', expected(:, i)');
            p_values(i) = p;
            chi2_stats{i} = stats;
        else
            p_values(i) = NaN;  % Not enough data to test
        end
    end

    rowLabels = {'Perivascular', 'Non-perivascular'};
    colLabels = ["↑↑", "↑→", "↑↓", "→↑", "→→", "→↓", "↓↑", "↓→", "↓↓"];

    residuals = (directionCounts - expected) ./ sqrt(expected);

%     % Plot mosaic-like figure
%     figure;
%     hold on;
%     total = sum(directionCounts(:));
%     x = 0;
% 
%     for j = 1:size(directionCounts, 2)
%         colTotal = sum(directionCounts(:, j));
%         width = colTotal / total;
%         y = 0;
% 
%         for i = 1:size(directionCounts, 1)
%             height = directionCounts(i, j) / colTotal;
%             color = [0.7 0.7 0.7]; % default gray
% 
%             if residuals(i, j) > 2
%                 color = [1 0.4 0.4]; % red
%             elseif residuals(i, j) < -2
%                 color = [0.4 0.4 1]; % blue
%             end
% 
%             rectangle('Position', [x y width height], 'FaceColor', color, 'EdgeColor', 'k');
%             y = y + height;
%         end
% 
%         x = x + width;
%     end
% 
%     xlim([0 1]);
%     ylim([0 1]);
%     axis off;
%     title('Mosaic Plot of Clusters by Cell Type');
% 
%     % Add legend
%     patch(NaN, NaN, [1 0.4 0.4], 'DisplayName', 'Residual > 2');
%     patch(NaN, NaN, [0.4 0.4 1], 'DisplayName', 'Residual < -2');
%     patch(NaN, NaN, [0.7 0.7 0.7], 'DisplayName', 'Expected');
%     legend('Location', 'eastoutside');
end
