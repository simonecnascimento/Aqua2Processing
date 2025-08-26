% Read PNG
img = imread('V:\2photon\Simone\Simone_Macrophages\Pf4Ai162-11\230630_FOV3\AQuA2\Pf4Ai162-11_230630_FOV3_run1_reg_Z01_green_Substack(1-927)\risingMaps_CH1\35.png');

% Create figure without borders
figure('Color', 'w', 'Visible', 'off');
imshow(img);
axis off;
set(gca, 'Position', [0 0 1 1]); % Fill figure

% Save as PDF
print('-dpdf', '-painters', '-r300', 'Pf4Ai162-11_230630_FOV3_Cell1_35.pdf');

close;
