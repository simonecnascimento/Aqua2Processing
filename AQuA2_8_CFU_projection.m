
clear all;

% Set the directory for the experiment you need
AQuA2_res_dir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\res';
AQuA2_cfu_dir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\cfu';
saveDir = 'D:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\cfu_overlay_images';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Get all .mat files in the directory
FilesAll = dir(fullfile(AQuA2_res_dir, '*_AQuA2.mat')); 

for x = 1:numel(FilesAll)
    resFileName = FilesAll(x).name;
    resFilePath = fullfile(FilesAll(x).folder, resFileName);
    session1 = load(resFilePath);
    session1.datPro = rescale(mean(single(session1.res.datOrg1), 4));
    
    baseName = erase(resFileName, '_AQuA2.mat');
    cfuFileName = [baseName '_AQuA_res_cfu.mat'];
    cfuFilePath = fullfile(AQuA2_cfu_dir, cfuFileName);
    cfuStruct = load(cfuFilePath);
    session1.cfuInfo1 = cfuStruct.cfuInfo1;

    session2 = session1;

    overlap = 0.4;  % how large IoU be considered as the same CFU
    [H, W] = size(session1.datPro);
    cfu1 = session1.cfuInfo1;
    cfu2 = session2.cfuInfo1;
    nCFU1 = size(cfu1, 1);
    nCFU2 = size(cfu2, 1);
    regionThr = 0.1;
    
    cfuMap1 = false(H, W, nCFU1);
    cfuMap2 = false(H, W, nCFU2);
    
    for i = 1:nCFU1
        cfuMap1(:, :, i) = cfu1{i, 3} > regionThr;
    end
    for i = 1:nCFU2
        cfuMap2(:, :, i) = cfu2{i, 3} > regionThr;
    end

    % Pair
    pairs1 = zeros(nCFU1, 2);
    pairs2 = zeros(nCFU2, 2);
        
    cfuMapCheck2 = reshape(cfuMap2, [], nCFU2);
    for i = 1:nCFU1
        id1 = i;
        pix = find(cfu1{id1, 3} > regionThr);
        candidates = find(sum(cfuMapCheck2(pix, :),1));
        IoUs = zeros(1, numel(candidates));
        for j = 1:numel(candidates)
            id2 = candidates(j);
            pix2 = find(cfu2{id2, 3} > regionThr);
            pixIn = intersect(pix, pix2);
            pixUnion = union(pix, pix2);
            IoUs(j) = numel(pixIn) / numel(pixUnion);
        end
        [IoU, id2] = max(IoUs);
        id2 = candidates(id2);
        if IoU > overlap
            pairs1(id1,:) = [id2, IoU];
            pairs2(id2,:) = [id1, IoU];
        end 
    end
    
    %
    common = cell(0,1);
    only1 = cell(0,1);
    only2 = cell(0,1);
    for i = 1:nCFU1
        if pairs1(i, 1) > 0
            id1 = i;
            id2 = pairs1(i);
            common{numel(common) + 1, 1} = (cfu1{id1, 3} + cfu2{id2, 3}) / 2;
        else
            only1{numel(only1) + 1, 1} = cfu1{i, 3};
        end
    end
    for i = 1:nCFU2
        if pairs2(i, 1) == 0
            only2{numel(only2) + 1, 1} = cfu2{i, 3};
        end
    end
    
    %
    datPro = rescale(session1.datPro);
    ov = cat(3, datPro, datPro, datPro);
    for i = 1:nCFU1
        x = randi(255,[1,3]);
        while (x(1)>0.8*255 && x(2)>0.8*255 && x(3)>0.8*255) || sum(x)<255
            x = randi(255,[1,3]);
        end
        ov(:, :, 1) = ov(:, :, 1) + 0.8 * x(1) / 255 * cfu1{i,3};
        ov(:, :, 2) = ov(:, :, 2) + 0.8 * x(2) / 255 * cfu1{i,3};
        ov(:, :, 3) = ov(:, :, 3) + 0.8 * x(3) / 255 * cfu1{i,3};
    end
    figure;
    imshow(ov)

    % Save as high-resolution TIFF (best for ImageJ)
    saveFileNameTIFF = [baseName '_cfuOverlay.tiff'];
    saveFilePathTIFF = fullfile(saveDir, saveFileNameTIFF);
    
    % Export the current figure to TIFF
    exportgraphics(gcf, saveFilePathTIFF, 'Resolution', 300);
    
    % Or, as PNG:
    saveFileNamePNG = [baseName '_cfuOverlay.png'];
    saveFilePathPNG = fullfile(saveDir, saveFileNamePNG);
    
    exportgraphics(gcf, saveFilePathPNG, 'Resolution', 300);
    
    % Close figure
    close;

    saveFileName = [baseName '_cfuOverlay.eps'];
    saveFilePath = fullfile(saveDir, saveFileName);

    print('-depsc2', '-painters', '-r300', saveFilePath);
    close;
end

%     % Construct base name without extension
%     [~, baseName, ~] = fileparts(resFileName);
% 
%     % Create the output file name
%     saveFileName = [baseName '_cfuOverlay.png'];
%     saveFilePath = fullfile(saveDir, saveFileName);

%% SAVE as VECTOR

clear; clc;

% === USER PATHS ===
AQuA2_res_dir = 'V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\res';
AQuA2_cfu_dir = 'V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\cfu';
saveDir = 'V:\2photon\Simone\Simone_Macrophages\AQuA2_Results\fullCraniotomy\baseline\cfu_overlay_images_vector';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% === PROCESS ALL FILES ===
FilesAll = dir(fullfile(AQuA2_res_dir, '*_AQuA2.mat'));

for x = 1:numel(FilesAll)
    % --- Load projection ---
    resFileName = FilesAll(x).name;
    resFilePath = fullfile(FilesAll(x).folder, resFileName);
    session = load(resFilePath);
    datPro = rescale(mean(single(session.res.datOrg1), 4));  % rescaled 0-1

    % --- Load CFU ---
    baseName = erase(resFileName, '_AQuA2.mat');
    cfuFileName = [baseName '_AQuA_res_cfu.mat'];
    cfuFilePath = fullfile(AQuA2_cfu_dir, cfuFileName);
    cfuStruct = load(cfuFilePath);
    cfuInfo = cfuStruct.cfuInfo1;
    nCFU = size(cfuInfo, 1);

    % === Visualization ===
    figure('Color','w','Units','pixels','Position',[100 100 800 600]);
    hold on;
    axis equal;
    axis off;
    set(gca,'YDir','reverse');

    % Display grayscale projection
    imagesc(datPro);
    colormap(gray);
    alpha(0.8);

    % Overlay CFU contours
    regionThr = 0.1;
    rng(1); % for reproducibility
    for i = 1:nCFU
        bw = cfuInfo{i,3} > regionThr;
        B = bwboundaries(bw);
        colorRand = rand(1,3)*0.8 + 0.2; % avoid too dark/light
        for k = 1:numel(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'Color', colorRand, 'LineWidth', 1);
        end
    end

    title(strrep(baseName,'_','\_'),'FontSize',14);

    % === Save as VECTOR PDF ===
    saveFileNamePDF = [baseName '_cfuOverlay_vector.pdf'];
    saveFilePathPDF = fullfile(saveDir, saveFileNamePDF);
    print('-dpdf','-painters','-r300', saveFilePathPDF);

    % === Save as VECTOR EPS ===
    saveFileNameEPS = [baseName '_cfuOverlay_vector.eps'];
    saveFilePathEPS = fullfile(saveDir, saveFileNameEPS);
    print('-depsc2','-painters','-r300', saveFilePathEPS);

    % === Save as high-res TIFF for ImageJ ===
    saveFileNameTIFF = [baseName '_cfuOverlay_highres.tiff'];
    saveFilePathTIFF = fullfile(saveDir, saveFileNameTIFF);
    exportgraphics(gca, saveFilePathTIFF, 'Resolution', 300);

    close;
end

disp('✅ Done: Vector overlays saved as PDF, EPS, and TIFF for all files.');
