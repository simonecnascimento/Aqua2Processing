function [categories, classifiedTable] = classifyAllCategories(dataTable, type)
    % Classify rows into rising, decreasing, and noChange categories
    % Ensure there are at least two columns
    if size(dataTable, 2) < 2
        error('Table does not have enough columns for classification.');
    end

    if type == 1 
        % Logical conditions for classification
        isRising = dataTable(:, 2) > 1.25 * dataTable(:, 1); %1.25
        isDecreasing = dataTable(:, 2) < 0.75 * dataTable(:, 1);  %0.75
        isNoChange = ~(isRising | isDecreasing); 
    
        % Extract rows based on classification
        categories.increase = dataTable(isRising, :);
        categories.decrease = dataTable(isDecreasing, :);
        categories.noChange = dataTable(isNoChange, :);
    
        % Ensure that each category is a table (if not already)
        categories.increase = array2table(categories.increase);
        categories.decrease = array2table(categories.decrease);
        categories.noChange = array2table(categories.noChange);
    
        % Add the classification label in the 5th column for each category
        categories.increase.Classification = repmat({'increase'}, size(categories.increase, 1), 1);
        categories.decrease.Classification = repmat({'decrease'}, size(categories.decrease, 1), 1);
        categories.noChange.Classification = repmat({'noChange'}, size(categories.noChange, 1), 1);
    
    %     % Define the variable names (you can adjust these based on your data)
    %     variableNames = {'preCSD', 'duringCSD', 'postCSD', 'cellType', 'Classification'};
    %     
    %     % Set the variable names for each category
    %     categories.rising.Properties.VariableNames = variableNames;
    %     categories.decreasing.Properties.VariableNames = variableNames;
    %     categories.noChange.Properties.VariableNames = variableNames;
    
        % Combine the categories back into a single table
        classifiedTable = [categories.increase; categories.decrease; categories.noChange];

    else

        % PRE-TO-DURING Logical conditions for classification
        increasingAcute = dataTable(:, 2) > 2 * dataTable(:, 1);
        decreasingAcute = dataTable(:, 2) < 0.5 * dataTable(:, 1);
        noChangeAcute = ~(increasingAcute | decreasingAcute);
    
        % PRE-TO-POST Logical conditions for classification
        increasingChronic = dataTable(:, 3) > 2 * dataTable(:, 1);
        decreasingChronic = dataTable(:, 3) < 0.5 * dataTable(:, 1);
        noChangeChronic = ~(increasingChronic | decreasingChronic);

        % Extract rows based on ACUTE classification
        categories.acute.increase = array2table(dataTable(increasingAcute, :));
        categories.acute.decrease = array2table(dataTable(decreasingAcute, :));
        categories.acute.noChange = array2table(dataTable(noChangeAcute, :));

        % Extract rows based on CHRONIC classification
        categories.chronic.increase = array2table(dataTable(increasingChronic, :));
        categories.chronic.decrease = array2table(dataTable(decreasingChronic, :));
        categories.chronic.noChange = array2table(dataTable(noChangeChronic, :));

        % Add the classification label in the 5th column for each category
        categories.acute.increase.Classification = repmat({'increase'}, size(categories.acute.increase, 1), 1);
        categories.acute.decrease.Classification = repmat({'decrease'}, size(categories.acute.decrease, 1), 1);
        categories.acute.noChange.Classification = repmat({'noChange'}, size(categories.acute.noChange, 1), 1);

        categories.chronic.increase.Classification = repmat({'increase'}, size(categories.chronic.increase, 1), 1);
        categories.chronic.decrease.Classification = repmat({'decrease'}, size(categories.chronic.decrease, 1), 1);
        categories.chronic.noChange.Classification = repmat({'noChange'}, size(categories.chronic.noChange, 1), 1);
    
        % Combine the categories back into a single table
        classifiedTable.Acute = [categories.acute.increase; categories.acute.decrease; categories.acute.noChange];
        classifiedTable.Chronic = [categories.chronic.increase; categories.chronic.decrease; categories.chronic.noChange];
    
    end
end
