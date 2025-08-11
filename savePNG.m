% Specify the folder containing the .fig files
folderPath = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\figures\all cells (except multinucleated)\network_digraph\cell_digraphs';  % Replace with your folder path

% Get a list of all .fig files in the folder
figFiles = dir(fullfile(folderPath, '*.fig'));

% Loop through each .fig file
for i = 1:length(figFiles)
    % Construct the full file name of the .fig file
    figFilePath = fullfile(folderPath, figFiles(i).name);
    
    % Open the .fig file
    fig = openfig(figFilePath, 'new', 'visible');
    
    % Construct the .png file name (same name as the .fig but with .png extension)
    [~, fileName, ~] = fileparts(figFiles(i).name);
    pngFilePath = fullfile(folderPath, [fileName, '.png']);
    
    % Save the figure as a .png file
    saveas(fig, pngFilePath, 'png');
    
    % Close the figure
    close(fig);
end
