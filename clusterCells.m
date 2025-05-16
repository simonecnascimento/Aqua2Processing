
cluster1 = combinedTable_clusters(clusterID == 1, :); % x=1,2
cluster2 = combinedTable_clusters(clusterID == 2, :); % x=3 (Cell2),6,10,19,23,28,35
cluster3 = combinedTable_clusters(clusterID == 3, :); % x=1,2
cluster4 = combinedTable_clusters(clusterID == 4, :); % x=4 (Cell7), 17,19,26,27,30,32
cluster5 = combinedTable_clusters(clusterID == 5, :); 
cluster6 = combinedTable_clusters(clusterID == 6, :); % x=1,3,4
cluster7 = combinedTable_clusters(clusterID == 7, :); % x=3,7
cluster8 = combinedTable_clusters(clusterID == 8, :); 
cluster9 = combinedTable_clusters(clusterID == 9, :); % x=40 (Cell15),17,18,20:23,66,68,69,71,80,81

clusterTarget = cluster9;

% no shadow
for x = 1:size(clusterTarget,1)
    fig = figure;
    plot(clusterTarget.dFF{x,1})
    cellID = string(clusterTarget.("Cell ID")(x));
    name = extractBefore(clusterTarget.("fileNameColumn")(x), "_reg");
    title(strcat(cellID, " ", name));
    ylim([-0.5, 5]);
    xline(1854, '--r', 'LineWidth', 1.5);
    %saveas(fig, 'Cell15_Cluster9_Curve', 'epsc');
end

for x = 35:size(clusterTarget,1)
    fig = figure;
    hold on;

    % Get the trace and convert time to minutes
    trace = clusterTarget.dFF{x,1};
    time_sec = (0:length(trace)-1) * (900 / 927);  % sec
    time_min = time_sec / 60;                     % min

    % Set y-limits
    ylow = -0.5;
    yhigh = 5;
    ylim([ylow, yhigh]);
    xlim([0,61]);
    set(gca, 'TickDir', 'out');  % Y-axis ticks outside

    % Convert frame range 1855–1917 to time (in minutes)
    t1 = 1855 * (900/927) / 60;
    t2 = 1917 * (900/927) / 60;

    % Draw gray shaded region
    xShade = [t1 t2 t2 t1];
    yShade = [ylow ylow yhigh yhigh];
    fill(xShade, yShade, [0.8 0.8 0.8], 'EdgeColor', 'none');

    % Plot trace
    plot(time_min, trace, 'k');

    % Labels
    %xlabel('Time (min)');
    %ylabel('dF/F');
    
    % Title
    cellID = string(clusterTarget.("Cell ID")(x));
    name = extractBefore(clusterTarget.("fileNameColumn")(x), "_reg");
    title(strcat(cellID, " ", name));

    % Save as EPS
    saveas(fig, sprintf('Cell%s_Curve_shadow_60min', cellID), 'epsc');
    close(fig); % optional: close figure to save memory
end



for x = 35:size(clusterTarget,1)
    fig = figure;
    hold on;
    time_sec = (0:length(trace)-1) * (900 / 927);  % or 0.970 for clarity
    time_min = time_sec / 60;
    % Set fixed y-limits
    ylow = -0.5;
    yhigh = 5;
    ylim([ylow, yhigh]);
    set(gca, 'TickDir', 'out');  % Y-axis ticks outside
    
    % Add light gray shaded region from timepoint 1855 to 1917
    xShade = [1855, 1917, 1917, 1855];
    yShade = [ylow, ylow, yhigh, yhigh];
    fill(xShade, yShade, [0.85 0.85 0.85], 'EdgeColor', 'none');  % Light gray
    
    % Plot the calcium trace on top
    plot(clusterTarget.dFF{x,1}, 'k')
    
    % Add title
    cellID = string(clusterTarget.("Cell ID")(x));
    name = extractBefore(clusterTarget.("fileNameColumn")(x), "_reg");
    title(strcat(cellID, " ", name));
    
    % Save as EPS
    saveas(fig, 'Cell7_Cluster4_Curve_shadow2', 'epsc');
end