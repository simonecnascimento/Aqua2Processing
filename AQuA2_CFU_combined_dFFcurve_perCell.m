%% Get dFF curve per cell across runs

%load '_anaysis.mat' from experiment you want
clear all;
load('Pf4Ai162-1_221122_FOV1_reg_Z01_green_Substack (200-5599)_analysis.mat');
% Extract the i-th row of the column
targetCellValue = strtrim('Cell 2'); % Your target value in the 'Cell ID' column
row = find(strcmp(resultsFinal{:, 1}, targetCellValue)); % Find the row number for 'Cell ID' value
dFF_Data = resultsFinal.dFF{row};
save(strcat(fileTemp, '_dFFcurve_', targetCellValue), 'dFF_Data');

%% create a combined dFF data
clear all;
dFF_Data1 = load("Pf4Ai162-1_221122_FOV2_run1_reg_Z01_green_Substack(1-900)_dFFcurve_cell 8.mat");
dFF_Data2 = load("Pf4Ai162-1_221122_FOV2_run1_reg_Z01_green_Substack(901-1800)_dFFcurve_cell 2.mat");
combined_dFF_Data = [dFF_Data1.dFF_Data,dFF_Data2.dFF_Data];

% Plot the curve for the i-th row
combined_dFF = figure;
plot(combined_dFF_Data); % Plot the curve for the i-th row
xlabel('Frames');
ylabel('dFF');
xlim=900;
title('Curve Plot');

% imageData = imread('Pf4Ai162-1_221122_FOV2_run1_reg_Z01_green_Substack(1-900)_combined_cells_figure.png');
% 
% % Display the image on the same figure
% hold on;  % This keeps the existing plot and adds the image
% image('CData', imageData);
% %, 'XData', [min(x) max(x)], 'YData', [min(y) max(y)]);
% 
% % Add a legend if needed
% legend('Sin(x)', 'Image');
% 
% % Turn off hold to revert to normal plotting behavior
% hold off;

saveas(combined_dFF, 'combined_dFF_cells(8,2)');

