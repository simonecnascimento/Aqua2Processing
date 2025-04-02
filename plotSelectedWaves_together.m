function plotSelectedWaves_together(dataTable, toplot, savePath)
    % Define colors
    pinkC3 = [255, 128, 128] / 255;
    blueC8 = [85, 160, 251] / 255;
    
    % Set figure layout
    Nrow = 4;
    Ncol = 19; 
    spGrid = reshape(1:Nrow*Ncol, Ncol, Nrow)';
    xLabelStr = 'Time (sec)';
    yLabelStr = '\DeltaF/F0';
    leftOpt = {[0.05,0.04], [0.08,0.03], [0.08,0.02]};
    
    % Create figure
    GLMresultFig = figure('WindowState', 'normal', 'Color', 'w');
    
    % Loop through selected plots
    for x = 1:length(toplot)
        y = toplot(x);
        sp(x) = subtightplot(Nrow, Ncol, spGrid(x,1:Ncol-1), leftOpt{:});
        hold on;
        
        % Choose color
        if x == 1 || x == 2
            plot(dataTable.dFF{y,1}, 'Color', pinkC3, 'LineWidth', 0.2);
        else
            plot(dataTable.dFF{y,1}, 'Color', blueC8, 'LineWidth', 0.2);
        end
        
        xlim([0, 900]);
        ylim([-0.3, inf]);
        set(gca, 'TickDir', 'out', 'TickLength', [0.008,0], 'Box', 'off');
        
        % Label axes
        if x == length(toplot)
            xlabel(xLabelStr);
            ylabel(yLabelStr);
        else
            set(gca, 'XTickLabel', []);
        end
    end
    
    % Save the figure
    if nargin > 2 && ~isempty(savePath)
        saveas(GLMresultFig, savePath, 'epsc');
    end
end
