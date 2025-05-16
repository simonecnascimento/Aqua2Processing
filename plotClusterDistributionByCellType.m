function [clustersPercentage, roundedData] = plotClusterDistributionByCellType(clusters, directions, cells, savePath)
    % Convert counts from strings to numbers
    perivascularCounts = str2double(clusters(2, 1:end-1));
    perivascularTotal = sum(perivascularCounts);
    perivascularPercentage = perivascularCounts / perivascularTotal;

    nonPerivascularCounts = str2double(clusters(3, 1:end-1));
    nonPerivascularTotal = sum(nonPerivascularCounts);
    nonPerivascularPercentage = nonPerivascularCounts / nonPerivascularTotal;

    % percentage table with labels and percentages
    clustersPercentage = [directions; num2cell(perivascularPercentage); num2cell(nonPerivascularPercentage)];
    clustersPercentage = [clustersPercentage, cells];

    % matrix of percentages (rows = [perivascular; non-perivascular])
    clustersPercentage2 = [perivascularPercentage; nonPerivascularPercentage] * 100;

    % Round to the nearest integer
    roundedData = round(clustersPercentage2);

    % Ensure each row sums to 100
    for i = 1:size(roundedData, 1)
        rowSum = sum(roundedData(i, :));
        diff = 100 - rowSum;
        if diff ~= 0
            [~, idx] = max(abs(roundedData(i, :)));
            roundedData(i, idx) = roundedData(i, idx) + diff;
        end
    end

    % Create figure
    fig = figure;
    % Perivascular pie chart
    subplot(1, 2, 1);
    pie(perivascularPercentage, directions);
    title('Perivascular Cells');

    % Non-perivascular pie chart
    subplot(1, 2, 2);
    pie(nonPerivascularPercentage, directions);
    title('Non-Perivascular Cells');

    % Save figure if savePath is provided
%     if nargin >= 4 && ~isempty(savePath)
%         if ~exist(savePath, 'dir')
%             mkdir(savePath);
%         end
%         filename = fullfile(savePath, 'clusterPie.png');
%         saveas(fig, filename);
%         fprintf('Figure saved to: %s\n', filename);
%     end


end
