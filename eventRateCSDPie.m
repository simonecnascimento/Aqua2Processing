% Extract cell types (column 13)
cellTypes = combinedTable_complete{:, 13};

% Get unique cell types
[uniqueTypes, ~, typeIdx] = unique(cellTypes);

% Labels in desired order
labels = ["↑↑", "↑→", "↑↓", ...
          "→↑", "→→", "→↓", ...
          "↓↑", "↓→", "↓↓"];

% Loop through each cell type
for i = 1:numel(uniqueTypes)
    currType = uniqueTypes(i);
    
    % Group indices for this cell type
    groupForType = groupIndex(typeIdx == i);

    % Count group frequencies
    groupCounts = histcounts(groupForType, 1:10);
    
    % Remove zero entries for clean display
    nonZeroIdx = groupCounts > 0;
    displayCounts = groupCounts(nonZeroIdx);
    displayLabels = labels(nonZeroIdx);

    % Create pie chart
    figure;
    p = pie(displayCounts);

    % Apply custom labels inside slices
    for j = 2:2:length(p)
        p(j).String = displayLabels(j/2);
        p(j).FontSize = 12;
        p(j).HorizontalAlignment = 'center';
        p(j).FontWeight = 'bold';
    end

    title(['Event Rate Change - Cell Type: ' + string(currType)]);
end
