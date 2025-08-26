% Observed counts (rows = P, NP; cols = PreCSD, DuringCSD, PostCSD)
AcuteIncrease = [7 7; 3 4]; %(rows = P, NP; cols = Persistent Increase, Persistent Decrease)
PersistentIncrease = [7 11; 3 34];  %(rows = P, NP; cols = Acute Increase, Persistent No Change)
PersistentDecrease = [7 25; 4 110];  %(rows = P, NP; cols = Acute Increase, Persistent No Change)

% Chi-square test
[chi2stat, p, df, expected] = deal([]);

% Calculate expected counts
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

% Display results
disp('Expected counts:');
disp(expected);
disp('Standardized residuals:');
disp(stdResiduals);
disp(['Chi-square = ' num2str(chi2stat) ', df = ' num2str(df) ', p = ' num2str(p)]);

% Labels
rowLabels = {'P','NP'};
colLabels = {'Increase','No Change','Decrease'};

% Create heatmap
figure;
h = heatmap(colLabels, rowLabels, stdResiduals, ...
    'Colormap', redbluecmap, ... % requires redbluecmap from File Exchange, or use built-in colormap
    'ColorLimits', [-2 2], ...   % emphasize >|2| values
    'CellLabelFormat','%.2f');

% Title and annotations
title(sprintf('Standardized residuals (Chi^2 = %.3f, p = %.4f)', chi2stat, p));
xlabel('Time point');
ylabel('Cell type');

% Optional: grid lines bold for emphasis
h.GridVisible = 'on';
