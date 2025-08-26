function [chi2stat, p, df, expected, stdResiduals] = plotChi2Residuals(O, rowLabels, colLabels)
% plotChi2Residuals - Compute chi-square test, standardized residuals, and plot heatmap
%
% Syntax:
%   [chi2stat, p, df, expected, stdResiduals] = plotChi2Residuals(O, rowLabels, colLabels)
%
% Inputs:
%   O          - Observed counts matrix (rows = groups, cols = categories)
%   rowLabels  - Cell array of row labels (e.g., {'P','NP'})
%   colLabels  - Cell array of column labels (e.g., {'Increase','No Change','Decrease'})
%
% Outputs:
%   chi2stat     - Chi-square statistic
%   p            - p-value
%   df           - Degrees of freedom
%   expected     - Expected counts under independence
%   stdResiduals - Standardized Pearson residuals

    % Validate inputs
    if nargin < 3
        error('You must provide observed matrix, row labels, and column labels.');
    end

    % Expected counts
    rowTotals = sum(O, 2);
    colTotals = sum(O, 1);
    grandTotal = sum(rowTotals);
    expected = (rowTotals * colTotals) / grandTotal;

    % Chi-square statistic
    chi2stat = sum(((O - expected).^2) ./ expected, 'all');

    % Degrees of freedom
    df = (size(O,1)-1) * (size(O,2)-1);

    % p-value
    p = 1 - chi2cdf(chi2stat, df);

    % Standardized residuals
    stdResiduals = (O - expected) ./ sqrt(expected);

    % Display in Command Window
    disp('Observed counts:'); disp(O);
    disp('Expected counts:'); disp(expected);
    disp('Standardized residuals:'); disp(stdResiduals);
    fprintf('Chi-square = %.4f, df = %d, p = %.6f\n', chi2stat, df, p);

    % Plot heatmap
    figure;
    try
        cmap = redbluecmap; % Requires redbluecmap from File Exchange
    catch
        cmap = parula; % fallback
    end
    h = heatmap(colLabels, rowLabels, stdResiduals, ...
        'Colormap', cmap, ...
        'ColorLimits', [-2 2], ...
        'CellLabelFormat', '%.2f');
    
    title(sprintf('Standardized residuals (Chi^2 = %.3f, p = %.4f)', chi2stat, p));
    xlabel('Category');
    ylabel('Group');
    h.GridVisible = 'on';
end
