for x = 1:size(combinedTable_clusters, 1)
    figure;
    plot(combinedTable_clusters.dFF{x,1})
end


postCSD30min = eventsByCell_all(:,4);
eventCounts_30min = cellfun(@(x) numel(x), postCSD30min);
eventHz_CSD30min = eventHz_byCell;
eventHz_CSD30min_postCSD = eventHz_CSD30min(:,4);

postCSD45min = eventsByCell_all(:,4);
eventCounts_45min = cellfun(@(x) numel(x), postCSD45min);
eventHz_CSD45min = eventHz_byCell;
eventHz_CSD45min_postCSD = eventHz_CSD45min(:,4);

eventHzComparison = [eventHz_CSD30min_postCSD; eventHz_CSD45min_postCSD];

