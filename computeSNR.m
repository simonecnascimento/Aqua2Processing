function snr_values = computeSNR(dFF_all, poly_order, frame_size)
    % COMPUTESNR Computes SNR for each row in dFF_all using FFT and Savitzky-Golay filter
    %
    %   snr_values = computeSNR(dFF_all, poly_order, frame_size)
    %
    %   Inputs:
    %       dFF_all    - A matrix where each row is a time-series signal (cells x time points)
    %       poly_order - Polynomial order for the Savitzky-Golay filter
    %       frame_size - Frame size for the Savitzky-Golay filter (should be odd)
    %
    %   Output:
    %       snr_values - A column vector of SNR values for each row (cell)
    
    [numCells, ~] = size(dFF_all);
    snr_values = nan(numCells, 1);

    for i = 1:numCells
        signal = dFF_all(i, :);
        
        % Compute FFT (optional visualization step, can be removed)
        fft_signal = fft(signal);
%         freq = linspace(0, 1, length(signal));
% 
%         % Visualize FFT
%         figure; 
%         plot(freq, abs(fft_signal)); 
%         title(['FFT of Cell ' num2str(i)]);

        % Apply Savitzky-Golay filter
        filtered_signal = sgolayfilt(signal, poly_order, frame_size(i,1));

        % Compute noise as the difference between original and filtered signal
        noise = signal - filtered_signal;

        % Calculate power of signal and noise
        signal_power = rms(filtered_signal .^ 2);
        noise_power = rms(noise .^ 2);

        % Compute SNR in dB
        if noise_power > 0
            snr_values(i) = 10 * log10(signal_power / noise_power);
        else
            snr_values(i) = NaN; % Avoid division by zero
        end
    end
end

