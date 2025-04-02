function plot_selected_dFF_waves(combinedTable_NM, saveFolder, saveFlag, selectedCells)
    % Check if save folder exists, if not, create it
    if saveFlag && ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

    % Loop through the selected rows
    for i = 1:length(selectedCells)
        x = selectedCells(i); % Get the row index

        % Create a figure
        dFFwave = figure('Visible', 'off'); % Set to 'on' if you want to display figures

        % Plot the data
        plot(combinedTable_NM.dFF{x,1});
        hold on;

        % Extract relevant data
        locationValue = combinedTable_NM.("Cell location (0,perivascular;1,adjacent;2,none)")(x,1);
        numberOfEvents = combinedTable_NM.("Number of Events")(x,1);

        % Assign location label
        if locationValue == 0
            locationLabel = "Perivascular";
        elseif locationValue == 2
            locationLabel = "Non-Perivascular";
        else
            locationLabel = "Other";
        end

        % Set labels and figure properties
        ylabel(sprintf('Cell: %d, %s, Events: %d', x, locationLabel, numberOfEvents));
        box off;
        xlim([0, 900]);

        % Remove y-axis (optional)
        %set(gca, 'YColor', 'none', 'YTick', [], 'XColor', 'none', 'XTick', []);

        % Save figure if saveFlag is true
        if saveFlag
            filename = sprintf('dFFwave_%d.png', x);
            filePath = fullfile(saveFolder, filename);
            saveas(dFFwave, filePath);
            close(dFFwave); % Close after saving to free memory
        else
            set(dFFwave, 'Visible', 'on'); % Show figure if not saving
        end
    end
end
