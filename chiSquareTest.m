function chiSquareTest(data)
    % Perform a chi-square test on the provided data
    % Input: 
    %   - data: A matrix with two columns representing the two categorical variables
    %            Column 1: Index, Column 2: Location
    
    % Generate a contingency table
    [observed, row_labels, col_labels] = crosstab(data(:,1), data(:,2));

    % Calculate expected frequencies assuming independence
    n = sum(observed(:));  % Total count
    row_sums = sum(observed, 2);  % Sum of rows
    col_sums = sum(observed, 1);  % Sum of columns
    expected = (row_sums * col_sums) / n;  % Expected frequencies under the null hypothesis
    
    % Calculate the chi-square statistic
    chi2stat = sum((observed - expected).^2 ./ expected, 'all');
    
    % Degrees of freedom
    df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
    
    % Calculate p-value from chi-square distribution
    p_value = 1 - chi2cdf(chi2stat, df);
    
    % Display results
    disp(['Chi-Square Statistic: ', num2str(chi2stat)]);
    disp(['Degrees of Freedom: ', num2str(df)]);
    disp(['p-value: ', num2str(p_value)]);
    
    % Interpretation
    if p_value < 0.05
        disp('There is a significant association between Index and Location.');
    else
        disp('There is no significant association between Index and Location.');
    end

    % Plot observed vs expected frequencies as a heatmap
    figure;
    subplot(1,2,1); % Left plot for observed
    heatmap(row_labels, col_labels, observed, 'ColorbarVisible', 'on');
    title('Observed Frequencies');
    
    subplot(1,2,2); % Right plot for expected
    heatmap(row_labels, col_labels, expected, 'ColorbarVisible', 'on');
    title('Expected Frequencies');
    
    % Optionally, you can also plot the residuals (difference between observed and expected)
    residuals = (observed - expected) ./ sqrt(expected); % Standardized residuals
    figure;
    heatmap(row_labels, col_labels, residuals, 'ColorbarVisible', 'on');
    title('Standardized Residuals');
end
