function groupProbs = compute_Clusterwise_probabilities(data)
    % Computes and displays the probability of each column within defined groups
    % data: a matrix (e.g., 2x9) of counts
    
    % Define the groups of columns using numeric indices
    groups = {[1 2 3], [4 5 6], [7 8 9]};  % Indices for each group
    groupLabels = {"↑", "x", "↓"};  % Labels for first arrow
    secondArrows = {  % Second arrow labels for each group
        {"↑↑", "↑x", "↑↓"},   % Labels for the first group (↑)
        {"x↑", "xx", "x↓"},   % Labels for the second group (→)
        {"↓↑", "↓x", "↓↓"}    % Labels for the third group (↓)
    };
    
    rowLabels = {'Perivascular', 'Non-perivascular'};  % Labels for rows
    
    % Initialize the probability matrix
    groupProbs = zeros(size(data));
    
    % Calculate probabilities
    for g = 1:length(groups)
        cols = groups{g};  % Get column indices for the current group
        groupData = data(:, cols);  % Extract the data for the current group
        rowSums = sum(groupData, 2);  % Sum over each row
        rowSums(rowSums == 0) = 1;   % Avoid division by zero
        
        % Compute probabilities for each group
        groupProbs(:, cols) = groupData ./ rowSums;  % Calculate probabilities
    end
    
    % Display results
    for row = 1:size(data, 1)
        fprintf('\n%s:\n', rowLabels{row});
        
        % Display each group
        for g = 1:length(groups)
            fprintf('  %s group:\n', groupLabels{g});
            cols = groups{g};  % Get the column indices for the current group
            secondArrowGroup = secondArrows{g};  % Get the second arrows for the current group
            for i = 1:length(cols)
                colIdx = cols(i);  % Get the current column index
                % Get the probability for this column and the count
                prob = groupProbs(row, colIdx) * 100;
                count = data(row, colIdx);
                fprintf('    %s: %.2f%% (Count = %d)\n', secondArrowGroup{i}, prob, count);
            end
        end
    end


%     % Create PIE figures to display the probabilities for each group
%     for g = 1:length(groups)
%         cols = groups{g};  % Get the column indices for the current group
%         secondArrowGroup = secondArrows{g};  % Get the second arrows for the current group
%         
%         % Probability data for the group
%         probs = groupProbs(:, cols) * 100;  % Convert to percentage
%         
%         % Display Pie chart for each row (Perivascular and Non-perivascular)
%         for row = 1:size(data, 1)
%             figure;
%             % Plot the pie chart for the current group and row
%             pie(probs(row, :), secondArrowGroup);
%             title(sprintf('%s  %s Group', rowLabels{row}, groupLabels{g}));
%             % Display the percentage in the legend
%             legend(arrayfun(@(x) sprintf('%s: %.2f%%', secondArrowGroup{x}, probs(row, x)), 1:length(secondArrowGroup), 'UniformOutput', false));
%         end
%     end



    % Create figure for displaying grouped bar plots
%     for g = 1:length(groups)
%         cols = groups{g};  % Get the column indices for the current group
%         secondArrowGroup = secondArrows{g};  % Get the second arrows for the current group
%         
%         % Probability data for the group
%         probs = groupProbs(:, cols) * 100;  % Convert to percentage
%         
%         % Create a new figure for each group
%         figure;
%         
%         % Bar plot for each row (Perivascular and Non-perivascular)
%         hold on;
%         bar(probs', 'grouped');
%         set(gca, 'XTickLabel', secondArrowGroup);
%         
%         % Title and labels
%         title(sprintf('%s Group', groupLabels{g}));
%         xlabel('Second Arrow');
%         ylabel('Probability (%)');
%         legend({'Perivascular', 'Non-perivascular'}, 'Location', 'Best');
%         
%         % Display the probabilities above the bars
% %         for i = 1:size(probs, 2)
% %             for j = 1:size(probs, 1)
% %                 text(i, probs(j, i) + 2, sprintf('%.2f%%', probs(j, i)), 'HorizontalAlignment', 'center');
% %             end
% %         end
% 
%         % Perform Chi-square test for each group (columns in the group)
%         observed = groupData;  % Observed counts for the current group
%         expected = repmat(sum(observed, 2) / size(observed, 2), 1, size(observed, 2));  % Expected counts assuming equal distribution
%         
%         % Chi-square test for each row (Perivascular and Non-perivascular)
%         p_values = zeros(1, size(observed, 1));  % Initialize p-values for each row
%         for i = 1:size(observed, 1)
%             [~, p_values(i)] = chi2gof(observed(i, :), 'Expected', expected(i, :));
%         end
%         
%         % Display p-values above the bars
%         for i = 1:size(p_values, 2)
%             text(i, max(probs(:, i)) + 5, sprintf('p=%.3f', p_values(i)), 'HorizontalAlignment', 'center');
%         end
% 
%         
%         hold off;

    % Create a figure for stacked bar charts
    figure;
    
    % Plot stacked bar for each group (Perivascular vs Non-perivascular)
    for g = 1:length(groups)
        cols = groups{g};
        secondArrowGroup = secondArrows{g};  % Get second arrow labels for the current group
        
        % Extract probabilities for the current group
        perivascularProbs = groupProbs(1, cols);  % First row: Perivascular
        nonPerivascularProbs = groupProbs(2, cols);  % Second row: Non-perivascular
        
        % Normalize the probabilities to percentages (between 0 and 100)
        perivascularProbs = perivascularProbs * 100;
        nonPerivascularProbs = nonPerivascularProbs * 100;

        % Create subplot for each group
        subplot(2, 3, g);  % 2 rows, 3 columns of subplots
        
        % Create a stacked bar plot for the current group
        barData = [perivascularProbs; nonPerivascularProbs]';  % Data to stack
        
        % Create a stacked bar chart
        h = bar(barData', 'stacked');
        
        % Set X-axis labels as second arrow types
        set(gca, 'XTickLabel', {'P', 'NP'});
        
        % Set title and labels
        title(groupLabels{g});
        %xlabel('Second Arrow');
        ylabel('Percentage (%)');
        ylim([0,150])

        % Set the color for each segment of the bar based on second arrows
        h(1).FaceColor = [0.4 0.8 0.4];  % Perivascular
        h(2).FaceColor = [0.8 0.4 0.4];  % Non-perivascular

        % Chi-square test for each subgroup (second arrow in group)
        pText = '';
        for i = 1:length(cols)
            counts = [data(1, cols(i)), data(2, cols(i))];  % Raw counts
            expected = sum(counts) / 2 * [1 1];
            if all(counts > 0)
                [~, p] = chi2gof(1:2, 'Freq', counts, 'Expected', expected, 'Emin', 1);
                pText = [pText, sprintf('%s: p = %.3f\n', secondArrowGroup{i}, p)];
            else
                pText = [pText, sprintf('%s: n/a\n', secondArrowGroup{i})];
            end
        end
    
        % Add p-values text box outside the plot
        annotation('textbox', ...
            [0.085 + mod(g-1,3)*0.31, 0.75 - floor((g-1)/3)*0.45, 0.12, 0.15], ...
            'String', pText, 'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10);

        % Legend for the stacked bars
        legend(groupLabels, 'Location', 'Northeast');

        
    end
end
    