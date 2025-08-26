% Specify the folder containing the .fig files
folderPath = 'V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\figures\all cells (except multinucleated)\network_digraph\';

folderPath = 'V:\2photon\Simone\Simone_Macrophages\Pf4Ai162-12\230717_FOV5\AQuA2\Pf4Ai162-12_230717_FOV5_run1_reg_Z01_green_Substack(1-927)\risingMaps_CH1';

% Get a list of all .fig files in the folder
figFiles = dir(fullfile(folderPath, '*.fig')); %.fig

% Loop through each .fig file
for i = 1:length(figFiles)
    % Construct the full file name of the .fig file
    figFilePath = fullfile(folderPath, figFiles(i).name);
    
    % Open the .fig file invisibly
    fig = openfig(figFilePath, 'invisible'); 

    % Optional: set paper size to match figure size for Illustrator
    set(fig, 'PaperPositionMode', 'auto');

    % Construct the .pdf file name
    [~, fileName, ~] = fileparts(figFiles(i).name);
    pdfFilePath = fullfile(folderPath, [fileName, '2.pdf']);

    % Save as PDF (vector format, best for Illustrator)
    print(fig, pdfFilePath, '-dpdf', '-painters');

    % Optional: Save as SVG if you prefer
    % svgFilePath = fullfile(folderPath, [fileName, '.svg']);
    % print(fig, svgFilePath, '-dsvg', '-painters');

    % Close the figure
    close(fig);
end
