function [best_filtered_signal, best_params] = optimizeSgolay(signal, Fs)
    % Function to optimize Savitzky-Golay filtering based on FFT and SNR analysis
    % Inputs:
    % - signal: Raw time series data
    % - Fs: Sampling frequency
    % Outputs:
    % - best_filtered_signal: Filtered signal with optimal parameters
    % - best_params: Struct containing best polynomial order and frame size
    
    % Step 1: Plot FFT of raw signal
    plotFFT(signal, Fs);
    
    % Step 2: Define candidate filter parameters
    frame_sizes = [11, 21, 41, 61];  % Ensure these are odd numbers
    poly_orders = [2, 3, 4];  % Test quadratic, cubic, and quartic fits
    
    best_SNR = -Inf;
    best_filtered_signal = signal;
    best_params = struct('poly_order', 0, 'frame_size', 0);
    
    % Step 3: Try different combinations and compute SNR
    for p = poly_orders
        for f = frame_sizes
            if f >= p || mod(f, 2) == 0
                continue; % Frame size must be greater than poly order
            end
            
            % Apply Savitzky-Golay filter
            filtered_signal = sgolayfilt(signal, p, f);
            
            % Compute SNR
            SNR_value = computeSNR(signal, filtered_signal);
            
            % Store the best filter settings
            if SNR_value > best_SNR
                best_SNR = SNR_value;
                best_filtered_signal = filtered_signal;
                best_params.poly_order = p;
                best_params.frame_size = f;
            end
        end
    end
    
    % Step 4: Plot and compare results
    figure;
    subplot(2,1,1);
    plot(signal, 'k'); hold on;
    plot(best_filtered_signal, 'r', 'LineWidth', 1.5);
    title(sprintf('Optimized Filtering (Poly Order: %d, Frame Size: %d)', ...
        best_params.poly_order, best_params.frame_size));
    legend('Raw Signal', 'Filtered Signal');
    
    % Step 5: Plot FFT of the filtered signal
    subplot(2,1,2);
    plotFFT(best_filtered_signal, Fs);
end

% Function to compute and plot FFT
function plotFFT(signal, Fs)
    N = length(signal);
    f = (0:N-1)*(Fs/N);  % Frequency axis
    fft_signal = abs(fft(signal)); % Compute FFT magnitude

    figure;
    plot(f(1:N/2), fft_signal(1:N/2)); % Plot half of the spectrum
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('FFT of Signal');
    grid on;
end

% Function to compute SNR
function SNR_value = computeSNR(original, filtered)
    noise = original - filtered;  % Estimate noise
    signal_power = rms(filtered)^2;
    noise_power = rms(noise)^2;
    SNR_value = 10 * log10(signal_power / noise_power); % SNR in dB
end
