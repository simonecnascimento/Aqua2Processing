function [h, p, stats] = chi2gof_from_table(tbl)
    % Total sum
    N = sum(tbl(:));
    
    % Row and column sums
    rowSums = sum(tbl, 2);
    colSums = sum(tbl, 1);
    
    % Expected counts under independence
    expected = rowSums * colSums / N;

    % Chi-square statistic
    chi2stat = sum((tbl - expected).^2 ./ expected, 'all');

    % Degrees of freedom = (#rows - 1) * (#cols - 1)
    df = (size(tbl,1)-1) * (size(tbl,2)-1);

    % p-value
    p = 1 - chi2cdf(chi2stat, df);

    % Output
    h = p < 0.05;
    stats.chi2stat = chi2stat;
    stats.df = df;
    stats.expected = expected;
end
