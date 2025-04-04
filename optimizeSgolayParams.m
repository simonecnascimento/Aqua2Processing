% function [best_poly_order, best_frame_size, best_snr, best_filtered_signals] = optimizeSgolayParams(dFF_all, poly_orders, frame_sizes)
%     % Optimize polynomial order and frame size for Savitzky-Golay filter
%     %
%     % Inputs:
%     %   dFF_all     - Matrix (cells x time points) containing fluorescence signals
%     %   poly_orders - Array of polynomial orders to test (e.g., 2:5)
%     %   frame_sizes - Array of frame sizes to test (must be odd and > poly order)
%     %
%     % Outputs:
%     %   best_poly_order  - Optimal polynomial order
%     %   best_frame_size  - Optimal frame size
%     %   best_snr         - Best SNR achieved
%     %   best_filtered_signals - Filtered signals with the best parameters
% 
%     [numCells, numTimePoints] = size(dFF_all);
%     best_snr = -Inf * ones(numCells, 1);  
%     best_poly_order = zeros(numCells, 1);
%     best_frame_size = zeros(numCells, 1);
%     best_filtered_signals = zeros(size(dFF_all));
% 
%     % Iterate over each cell
%     for i = 1:numCells
%         signal = dFF_all(i, :);
% 
%         % Iterate through polynomial orders
%         for p = poly_orders
%             % Iterate through frame sizes (must be greater than poly order)
%             for f = frame_sizes
%                 if f <= p
%                     continue; % Skip invalid frame sizes
%                 end
% 
%                 % Apply Savitzky-Golay filter
%                 filtered_signal = sgolayfilt(signal, p, f);
%                 noise = signal - filtered_signal;
%                 snr_value = var(filtered_signal) / var(noise);
% 
%                 % Update best parameters if current SNR is higher
%                 if snr_value > best_snr(i)
%                     best_snr(i) = snr_value;
%                     best_poly_order(i) = p;
%                     best_frame_size(i) = f;
%                     best_filtered_signals(i, :) = filtered_signal;
%                 end
%             end
%         end
%     end
% end

%%
function [best_poly_order, best_frame_size, best_snr, best_filtered_signals] = optimizeSgolayParams(dFF_all, poly_orders, frame_sizes)
    % Optimize polynomial order and frame size for Savitzky-Golay filter
    %
    % Inputs:
    %   dFF_all     - Matrix (cells x time points) containing fluorescence signals
    %   poly_orders - Array of polynomial orders to test (e.g., 2:5)
    %   frame_sizes - Array of frame sizes to test (must be odd and > poly order)
    %
    % Outputs:
    %   best_poly_order  - Optimal polynomial order
    %   best_frame_size  - Optimal frame size
    %   best_snr         - Best SNR achieved
    %   best_filtered_signals - Filtered signals with the best parameters

    [numCells, numTimePoints] = size(dFF_all);
    best_snr = -Inf * ones(numCells, 1);  
    best_poly_order = zeros(numCells, 1);
    best_frame_size = zeros(numCells, 1);
    best_filtered_signals = zeros(size(dFF_all));

    % Parallelize the loop over cells
    parfor i = 1:numCells
        signal = dFF_all(i, :);
        best_snr_i = -Inf;  % Best SNR for the current signal
        best_poly_order_i = 0; % Best poly order for the current signal
        best_frame_size_i = 0; % Best frame size for the current signal
        best_filtered_signal_i = zeros(1, numTimePoints);  % To store the filtered signal

        % Iterate through polynomial orders
        for p = poly_orders
            % Iterate through frame sizes (must be greater than poly order)
            for f = frame_sizes
                if f <= p
                    continue; % Skip invalid frame sizes
                end

                % Apply Savitzky-Golay filter
                filtered_signal = sgolayfilt(signal, p, f);
                noise = signal - filtered_signal;
                snr_value = var(filtered_signal) / var(noise);

                % Update best parameters if current SNR is higher
                if snr_value > best_snr_i
                    best_snr_i = snr_value;
                    best_poly_order_i = p;
                    best_frame_size_i = f;
                    best_filtered_signal_i = filtered_signal;
                end
            end
        end

        % Assign the best values for this cell to the main arrays
        best_snr(i) = best_snr_i;
        best_poly_order(i) = best_poly_order_i;
        best_frame_size(i) = best_frame_size_i;
        best_filtered_signals(i, :) = best_filtered_signal_i;
    end
end

