function cleanDataDistances_um = cellPairsDistanceDistribution(allUpperTriValues)
    % Function to calculate and visualize the distribution of distances
    % between cells with concurrent activity.
    %
    % Inputs:
    %   allUpperTriValues - Cell array where each cell contains a vector of
    %     distance values for different experiments or conditions.
    %
    % Outputs:
    %   cleanDataDistances_um - Cleaned array of all distance values (NaNs removed).
    
    % Initialize an empty array to store all distance values
    allValuesDistances_um = [];
    
    % Loop through each cell in the input array
    for k = 1:length(allUpperTriValues)
        currentVector = allUpperTriValues{k};
        
        % Append the current vector's values to the master array
        allValuesDistances_um = [allValuesDistances_um; currentVector]; 
    end
    
    % Remove NaN values from the combined data
    cleanDataDistances_um = allValuesDistances_um(~isnan(allValuesDistances_um));
    cleanDataDistances_um_unique = unique(cleanDataDistances_um);   

    % Plot the histogram of the cleaned data
    figure;
    histogram(cleanDataDistances_um);
    title('Distance Metrics for Cells with Concurrent Activity');
    xlabel('Distance (um)');
    ylabel('Number of Cell Pairs');
    
    % Optionally return the cleaned data for further analysis
end
