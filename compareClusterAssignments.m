function [matchIdx, diffIdx, clusterTable, crosstabMatrix] = compareClusterAssignments(clusterCol1, clusterCol2)
%COMPARECLUSTERASSIGNMENTS Compare and associate two sets of cluster assignments.
%
% INPUTS:
%   clusterCol1 - Nx1 vector of cluster IDs (e.g., from method A)
%   clusterCol2 - Nx1 vector of cluster IDs (e.g., from method B)
%
% OUTPUTS:
%   matchIdx        - logical Nx1 array where clusterCol1 == clusterCol2
%   diffIdx         - indices where clusterCol1 ~= clusterCol2
%   clusterTable    - table showing the cluster assignment pairs
%   crosstabMatrix  - matrix showing how many times each pair of clusters co-occur

    cluster_eventRate = combinedTable_clusters.eventRate_clusterID;
    cluster_dFF = combinedTable_clusters.dFF_clusterID;

    % Check inputs
    if numel(cluster_eventRate) ~= numel(cluster_dFF)
        error('Inputs must be the same length.');
    end

    % Ensure column vectors
    cluster_eventRate = cluster_eventRate(:);
    cluster_dFF = cluster_dFF(:);

    % Match and mismatch indices
    matchIdx = cluster_eventRate == cluster_dFF;
    diffIdx = find(~matchIdx);

    % Create table of associations
    clusterTable = table(cluster_eventRate, cluster_dFF, matchIdx, ...
        'VariableNames', {'Cluster1', 'Cluster2', 'IsMatch'});

    % Crosstab (co-occurrence matrix)
    [~,~,group1] = unique(cluster_eventRate);
    [~,~,group2] = unique(cluster_dFF);
    crosstabMatrix = crosstab(group1, group2);

    % Optional: Display summary
    fprintf('Matching assignments: %d of %d (%.1f%%)\n', ...
        sum(matchIdx), numel(matchIdx), 100 * mean(matchIdx));

    % Optional: Plot heatmap of co-occurrence
    figure;
    heatmap(crosstabMatrix, ...
        'XLabel', 'Cluster dFF', 'YLabel', 'Cluster event rate', ...
        'Title', 'Cluster Co-Occurrence Heatmap');
end
