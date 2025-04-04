function tbl = chiSquaredTestForAssociation(data)
    % Create a contingency table between cell location and cluster index
    [tbl, clusterLabels, locLabels] = crosstab(data(:,1), data(:,2));
    
    % Display contingency table
    disp('Contingency Table:');
    disp(tbl);

    % Perform Chi-squared test for independence
    [chi2stat, p_value, df] = chi2gof(tbl(:), 'Ctrs', 1:numel(tbl), 'Expected', mean(tbl(:)) * ones(size(tbl(:))));

    % Display the results
    fprintf('Chi-squared statistic: %.4f\n', chi2stat);
    fprintf('p-value: %.4f\n', p_value);
    fprintf('Degrees of Freedom: %.0f\n', df.df);

    if p_value < 0.05
        disp('There is a significant association between cell location and cluster index.');
    else
        disp('There is no significant association between cell location and cluster index.');
    end

    % Plot the contingency table as a heatmap
    figure;
    imagesc(tbl);  % Create a heatmap
    colorbar;  % Show color bar
    title('Contingency Table Heatmap');
    xlabel('Cell Location');
    ylabel('Cluster');
    %xticks(1:length(locLabels));  % Set x-ticks to the cluster labels
    %xticklabels(locLabels);
    %yticks(1:length(clusterLabels));  % Set y-ticks to the location labels
    %yticklabels(clusterLabels);
    colormap parula;  % Choose a color map for the heatmap

     % Add numbers inside squares with white text and black background
    for i = 1:size(tbl, 1)  % Loop over rows (clusters)
        for j = 1:size(tbl, 2)  % Loop over columns (locations)
            text(j, i, num2str(tbl(i, j)), 'HorizontalAlignment', 'center', ...
                'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'k');
        end
    end

  
end
