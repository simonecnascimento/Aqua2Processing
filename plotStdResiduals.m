function plotStdResiduals(matrix, rowLabels, colLabels, outputDir)
% plotStdResiduals - Plots a heatmap of standardized Pearson residuals for a contingency table
%
% Inputs:
%   All_matrix - contingency table (matrix)
%   rowLabels  - cell array of row labels
%   colLabels  - cell array of column labels
% Example:
%   plotStdResiduals(All_matrix, {'Acute ↑', 'Acute x', 'Acute ↓'}, {'Chronic ↑', 'Chronic x', 'Chronic ↓'});

    % Calculate expected counts
    rowSums = sum(matrix, 2);   % Sum across rows
    colSums = sum(matrix, 1);   % Sum across columns
    grandTotal = sum(matrix(:));
    expected = (rowSums * colSums) / grandTotal;  % Outer product

    % Calculate standardized residuals
    stdResiduals = (matrix - expected) ./ sqrt(expected);

    % Create figure
    figure;
    h = heatmap(colLabels, rowLabels, stdResiduals, ...
        'Colormap', redbluecmap, ...
        'ColorLimits', [-max(abs(stdResiduals(:))), max(abs(stdResiduals(:)))]); % symmetric color scaling
    
    % Add title and labels
    h.XLabel = 'Chronic Response';
    h.YLabel = 'Acute Response';

    %save
    fullPath = fullfile(outputDir, 'Standard residuals.eps');
    print(gcf, fullPath, '-depsc', '-r300');
end


