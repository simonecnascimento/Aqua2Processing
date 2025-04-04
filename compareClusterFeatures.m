function compareClusterFeatures(clusterFeatures, clusterFeatures_cluster1, clusterFeatures_cluster2)
    % This function performs Mann-Whitney U tests between two clusters for each feature,
    % computes the median and IQR for each cluster, and creates bar plots with significance annotations.
    
    % Get the feature names (excluding 'Cluster' column)
    featureNames = clusterFeatures.Properties.VariableNames(2:end);
    numFeatures = length(featureNames);

    % Initialize p-values for the Mann-Whitney U tests
    p_values = zeros(numFeatures, 1);
    testResults = table(featureNames', p_values, 'VariableNames', {'Feature', 'P_Value'});

    % Perform Mann-Whitney U test (ranksum) for each feature
    for i = 1:numFeatures
        data1 = clusterFeatures_cluster1{:, featureNames{i}}; % Data for Cluster 1
        data2 = clusterFeatures_cluster2{:, featureNames{i}}; % Data for Cluster 2
        
        % Mann-Whitney U test
        p_values(i) = ranksum(data1, data2);  
    end

    % Store p-values in the table
    testResults.P_Value = p_values;

    % Display p-values for each feature
    disp(testResults);

    % Compute Median and IQR for each feature
    medianValues = zeros(numFeatures, 2);
    iqrValues = zeros(numFeatures, 2);

    for i = 1:numFeatures
        % Extract feature data
        data1 = clusterFeatures_cluster1{:, featureNames{i}}; % Cluster 1
        data2 = clusterFeatures_cluster2{:, featureNames{i}}; % Cluster 2
        
        % Compute Median and IQR
        medianValues(i, 1) = median(data1);
        medianValues(i, 2) = median(data2);
        
        iqrValues(i, 1) = iqr(data1);
        iqrValues(i, 2) = iqr(data2);
    end

    % Plot subplots with bar graphs
    figure;
    for i = 1:numFeatures
        subplot(ceil(numFeatures / 3), 3, i); % Arrange subplots in rows of 3
        hold on;
        
        % Bar plot for median values
        barHandle = bar(medianValues(i, :)); 

        % Set different colors for clusters
        barHandle.CData(1, :) = [1 0.5 0]; % Orange for Cluster 1
        barHandle.CData(2, :) = [0 0.6 0]; % Green for Cluster 2

        % Error bars (IQR)
        x = [1, 2]; % X positions of bars
        errorbar(x, medianValues(i, :), iqrValues(i, :), 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

        % Add p-value annotation
        if p_values(i) < 0.001
            sigText = '***';  % Highly significant
        elseif p_values(i) < 0.01
            sigText = '**';   % Very significant
        elseif p_values(i) < 0.05
            sigText = '*';    % Significant
        else
            sigText = 'n.s.'; % Not significant
        end
        
        % Display p-value above bars
        text(1.5, max(medianValues(i, :)) + max(iqrValues(i, :)) * 0.1, sigText, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

        % Formatting
        title(featureNames{i}, 'Interpreter', 'none');
        xticks([1 2]); 
        xticklabels({'Cluster 1', 'Cluster 2'});
        ylabel('Median ± IQR');
        grid on;
        hold off;
    end

    % Super title for all subplots
    sgtitle('Comparison of Features Between Clusters'); 
end
