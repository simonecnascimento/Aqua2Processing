function best_frame_sizes = optimizeFrameSize(dFF_all, poly_order, frame_range)
    % OPTIMIZEFRAMESIZE Finds the best frame size for each cell in dFF_all
    %
    %   best_frame_sizes = optimizeFrameSize(dFF_all, poly_order, frame_range)
    %
    %   Inputs:
    %       dFF_all    - Matrix where each row is a time-series signal (cells x time points)
    %       poly_order - Polynomial order for Savitzky-Golay filter
    %       frame_range - A vector of possible frame sizes (must be odd numbers)
    %
    %   Output:
    %       best_frame_sizes - A column vector of best frame sizes for each row (cell)

    [numCells, numTimePoints] = size(dFF_all);
    best_frame_sizes = nan(numCells, 1);  % Initialize
    best_filtered_signal = nan(size(dFF_all)); % Pre-allocate

    parfor i = 1:numCells
        signal = dFF_all(i, :);
        best_snr = -Inf;
        best_frame_size = frame_range(1);
        best_filtered = nan(1, size(dFF_all, 2)); % Preallocate per cell

        for f_size = frame_range
            if f_size > numTimePoints  % Skip invalid sizes
                continue;
            end

            % Apply Savitzky-Golay filter
            filtered_signal = sgolayfilt(signal, poly_order, f_size);
            noise = signal - filtered_signal;

            % Compute SNR
            noise_var = var(noise);
            if noise_var == 0
                continue;
            end
            snr_value = var(filtered_signal) / noise_var;

            % Select best frame size
            if snr_value > best_snr
                best_snr = snr_value;
                best_frame_size = f_size;
                best_filtered = filtered_signal; % Store the best signal
            end

%             % Stop early if improvement is < 1%
%             if abs(snr_value - prev_snr) / prev_snr < 0.01
%                 break;
%             end
%             prev_snr = snr_value;

        end

        best_frame_sizes(i) = best_frame_size;
        best_filtered_signal(i, :) = best_filtered; % Assign the best signal safely

        %plotFilteredSignal(dFF_all(i, :), best_filtered_signal(i, :), best_frame_sizes(i), numTimePoints, i);
    end
end
