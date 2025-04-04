function snr_values = computeSNR(dFF_all, method, param)
    % COMPUTESNR Computes SNR using different methods
    % 
    %   snr_values = computeSNR(dFF_all, method, param)
    % 
    %   Inputs:
    %       dFF_all - A matrix where each row is a time-series signal (cells x time points)
    %       method  - A string indicating the SNR calculation method:
    %                 'ma'  -> Moving Average filtering (param = window size)
    %                 'fft' -> FFT-based noise estimation (param = cutoff frequency ratio)
    %       param   - Numeric parameter based on method
    % 
    %   Output:
    %       snr_values - A column vector of SNR values for each row (cell)

    [numCells, numTimePoints] = size(dFF_all);
    snr_values = nan(numCells, 1);

    switch method
        case 'ma'  % Moving Average Method
            window_size = param;
            for i = 1:numCells
                signal = dFF_all(i, :);

                % Apply moving average filter
                smooth_signal = movmean(signal, window_size);

                % Compute noise as the difference between original and smooth signal
                noise = signal - smooth_signal;

                % Calculate power of signal and noise
                signal_power = mean(smooth_signal .^ 2);
                noise_power = mean(noise .^ 2);

                % Compute SNR in dB
                if noise_power > 0
                    snr_values(i) = 10 * log10(signal_power / noise_power);
                else
                    snr_values(i) = NaN;
                end
            end

        case 'fft'  % FFT-Based Method
            cutoff_freq_ratio = param;
            for i = 1:numCells
                signal = dFF_all(i, :);

                % Compute FFT
                fft_signal = fft(signal);
                power_spectrum = abs(fft_signal).^2;

                % Define cutoff index
                cutoff_idx = round(numTimePoints * cutoff_freq_ratio);

                % Estimate signal power from low frequencies
                signal_power = sum(power_spectrum(1:cutoff_idx));

                % Estimate noise power from high frequencies
                noise_power = sum(power_spectrum(cutoff_idx+1:end));

                % Compute SNR in dB
                if noise_power > 0
                    snr_values(i) = 10 * log10(signal_power / noise_power);
                else
                    snr_values(i) = NaN;
                end
            end

        otherwise
            error('Invalid method. Use "ma" for Moving Average or "fft" for FFT-based SNR calculation.');
    end
end




% function snr_values = computeSNR(dFF_all, poly_order, frame_size)
%     % COMPUTESNR Computes SNR for each row in dFF_all using FFT and Savitzky-Golay filter
%     %
%     %   snr_values = computeSNR(dFF_all, poly_order, frame_size)
%     %
%     %   Inputs:
%     %       dFF_all    - A matrix where each row is a time-series signal (cells x time points)
%     %       poly_order - Polynomial order for the Savitzky-Golay filter
%     %       frame_size - Frame size for the Savitzky-Golay filter (should be odd)
%     %
%     %   Output:
%     %       snr_values - A column vector of SNR values for each row (cell)
%     
%     [numCells, ~] = size(dFF_all);
%     snr_values = nan(numCells, 1);
% 
%     for i = 1:numCells
%         signal = dFF_all(i, :);
%         
%         % Compute FFT (optional visualization step, can be removed)
%         fft_signal = fft(signal);
% %         freq = linspace(0, 1, length(signal));
% % 
% %         % Visualize FFT
% %         figure; 
% %         plot(freq, abs(fft_signal)); 
% %         title(['FFT of Cell ' num2str(i)]);
% 
%         % Apply Savitzky-Golay filter
%         filtered_signal = sgolayfilt(signal, poly_order, frame_size(i,1));
% 
%         % Compute noise as the difference between original and filtered signal
%         noise = signal - filtered_signal;
% 
%         % Calculate power of signal and noise
%         signal_power = rms(filtered_signal .^ 2);
%         noise_power = rms(noise .^ 2);
% 
%         % Compute SNR in dB
%         if noise_power > 0
%             snr_values(i) = 10 * log10(signal_power / noise_power);
%         else
%             snr_values(i) = NaN; % Avoid division by zero
%         end
%     end
% end
% 
