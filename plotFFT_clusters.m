poly_order = 6;  % Range of polynomial orders for Savitzky-Golay filter
frame_size = 7;
darkPurple = [0.5, 0, 0.5];  % Dark purple
darkGreen = [0, 0.5, 0];     % Dark green
colors = [darkPurple; darkGreen];  % Colors for Cluster 1 and 2
figure;
for k = 1:2
    % 1. Detrend the signal (remove systematic drifts like tissue drift, bleaching)
    detrended_signal = detrend(dFF_all(idx == k, :));
    
    % 2. Apply Savitzky-Golay filter for smoothing (window size 5, polynomial order 3)
    smoothed_signal = sgolayfilt(detrended_signal, poly_order, frame_size);
    
    % 3. Fourier Transformation: Calculate the power spectrum density (PSD)
    L = length(smoothed_signal);  % Length of the signal
    
    % Compute FFT
    Y = fft(smoothed_signal);

    if k == 2
        % Low pass filter
        Fc = 0.1; % Cutoff frequency in Hz
        order = 6; % Filter order (higher order gives steeper roll-off)
        % Create the low-pass filter
        [b, a] = butter(order, Fc/(Fs/2), 'low'); 
        % Assuming x is your signal
        Y = filter(b, a, Y);
    end

    P2 = abs(Y/L);                
    P1 = P2(1:L/2+1);            
    P1(2:end-1) = 2*P1(2:end-1);  % Adjust for single-sided spectrum
    
    % Frequency axis
    f = Fs * (0:(L/2)) / L;
    
    % Compute Power Spectrum Density (PSD)
    psd = P1.^2;
    
    % Plot the frequency spectrum for high-peak data
    subplot(2, 1, k);
    plot(f, psd, 'Color', colors(k,:), LineWidth=1);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    %ylim([0, 0.003])
    xlim([-0.005,0.52]);
    %title('Frequency Spectrum of Cluster');
    grid off;
    box off;
end

%%
% Example: Plot FFT for high-peak data (n=463)
signal_high_peak = dFF_all(idx == 2, :);
Fs = 1.03;  % Sampling frequency (Hz)
L = length(signal_high_peak);
Y_high = fft(signal_high_peak);  % Compute FFT

% Low pass filter
% Define parameters
Fc = 0.1; % Cutoff frequency in Hz
order = 6; % Filter order (higher order gives steeper roll-off)
% Create the low-pass filter
[b, a] = butter(order, Fc/(Fs/2), 'low'); 
% Assuming x is your signal
Y_high_filtered = filter(b, a, Y_high);

P2_high = abs(Y_high_filtered/L).^2;  % Two-sided spectrum
P1_high = P2_high(1:L/2+1);  % Single-sided spectrum
P1_high(2:end-1) = 2*P1_high(2:end-1);  % Correct for single-sided spectrum

% Frequency vector
f = Fs*(0:(L/2))/L;

% Plot the frequency spectrum for high-peak data
figure;
subplot(2, 1, 1);
plot(f, psd);
xlabel('Frequency (Hz)');
ylabel('Power');
%ylim([0, 0.013])
xlim([0,0.52]);
title('Frequency Spectrum of Cluster 2');
grid on;

%% Example: Plot FFT for low-peak data (n=40)
signal_low_peak = dFF_all(idx == 1, :); 
Y_low = fft(signal_low_peak);  % Compute FFT

% % Low pass filter
% % Define parameters
% Fc = 0.1; % Cutoff frequency in Hz
% order = 6; % Filter order (higher order gives steeper roll-off)
% % Create the low-pass filter
% [b, a] = butter(order, Fc/(Fs/2), 'low'); 
% % Assuming x is your signal
% Y_low_filtered = filter(b, a, Y_low);


P2_low = abs(Y_low/L).^2;  % Two-sided spectrum
P1_low = P2_low(1:L/2+1);  % Single-sided spectrum
P1_low(2:end-1) = 2*P1_low(2:end-1);  % Correct for single-sided spectrum

% Plot the frequency spectrum for low-peak data
subplot(2, 1, 2);
plot(f, P1_low);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Frequency Spectrum of Cluster 1');
%ylim([0, 0.013])
xlim([0,0.52]);
grid on;

%% Low pass filter
% Define parameters
Fc = 0.5; % Cutoff frequency in Hz
order = 4; % Filter order (higher order gives steeper roll-off)
% Create the low-pass filter
[b, a] = butter(order, Fc/(Fs/2), 'low'); 
% Assuming x is your signal
x_filtered = filter(b, a, x);

%% Example: Plot FFT for all data
signalAll = dFF_all(:, :);  % Replace with actual data
Fs = 1.03;  % Sampling frequency (Hz)
L = length(signalAll);
Y_all = fft(signalAll);  % Compute FFT
P2_all = abs(Y_all/L).^2;  % Two-sided spectrum
P1_all = P2_all(1:L/2+1);  % Single-sided spectrum
P1_all(2:end-1) = 2*P1_all(2:end-1);  % Correct for single-sided spectrum

% Plot the frequency spectrum for low-peak data
figure;
plot(f, P1_all);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Frequency Spectrum of All Data');

%% Overlap frequency spectrum graphs

% Define custom colors
darkPurple = [0.5, 0, 0.5];  % Dark purple
darkGreen = [0, 0.5, 0];     % Dark green
colors = [darkPurple; darkGreen];  % Colors for Cluster 1 and 2

% Example: Plot FFT for high-peak data (n=463)
signal_high_peak = dFF_all(idx == 2, :);
Fs = 1.03;  % Sampling frequency (Hz)
L = length(signal_high_peak);
Y_high = fft(signal_high_peak);  % Compute FFT
P2_high = abs(Y_high/L).^2;  % Two-sided spectrum
P1_high = P2_high(1:L/2+1);  % Single-sided spectrum
P1_high(2:end-1) = 2*P1_high(2:end-1);  % Correct for single-sided spectrum

% Example: Plot FFT for low-peak data (n=40)
signal_low_peak = dFF_all(idx == 1, :); 
Y_low = fft(signal_low_peak);  % Compute FFT
P2_low = abs(Y_low/L).^2;  % Two-sided spectrum
P1_low = P2_low(1:L/2+1);  % Single-sided spectrum
P1_low(2:end-1) = 2*P1_low(2:end-1);  % Correct for single-sided spectrum

% Frequency vector
f = Fs*(0:(L/2))/L;

% Plot the frequency spectrum for both high-peak and low-peak data on the same plot
figure;
hold on; % Hold the current plot to overlay both spectra

% Plot the high-peak data (Cluster 2)
h2 = plot(f, P1_high, 'Color', colors(2,:), 'LineWidth', 2);  

% Plot the low-peak data (Cluster 1)
h1 = plot(f, P1_low, 'Color', colors(1,:), 'LineWidth', 2);  % Low peak in blue

% Set labels, title, and other properties
xlabel('Frequency (Hz)');
ylabel('Power');
title('Overlapping Frequency Spectra of High and Low Peak Data');
legend([h1, h2],{'Cluster 1', 'Cluster 2'}, 'Location', 'northeast');
ylim([0, 0.013]); % Adjust y-axis limits as needed
grid off; box off;
xlim([-0.005,0.52])
hold off;  % Release the hold to prevent further plotting on the same figure



